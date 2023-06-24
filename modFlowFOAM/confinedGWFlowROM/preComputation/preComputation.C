/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
/*-----------------------------------------------------------------------*\
Class
    Foam::preComputation

Group
    preComputation/confinedGWFlowROM/modFlowFOAM

Description
    Pre-computation required for reduced-order modeling of unsteady flow 
    through a homogeneous/heterogeneous confined aquifer
\*-----------------------------------------------------------------------*/
#include "fvCFD.H"
#include "fvOptions.H"
#include "IFstream.H"
#include "OFstream.H"
#include "readCSV.H"
#include "SVD.H"
#include "simpleMatrix.H"

int main(int argc, char *argv[])
{
    std::clock_t startT= std::clock(); 									//Start Time
    
    argList::addNote
    (
        "Inverse of the Projected Stiffness Matrix"
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "RectangularMatrix.H"
    #include "SquareMatrix.H"
    #include "scalarMatrices.H"
    
    //----------------------------------------------------------------------------------------------------//
	#include "readModelParameters.H" 									//Reading the input parameters
	
	if (conductivityFieldMode.match("T"))								//Transmissivity provided
	{
		T = T;
	}
	else if (conductivityFieldMode.match("HC"))							//Hydraulic Conductivity provided 
	{
		T = b*Ks;
	}
	
    //----------------------------------------------------------------------------------------------------//
	//Reading the snapshot matrix from CSV files into a Matrix
	RectangularMatrix <scalar> matSnap = readCSV("snapshotMatrix.csv");
	
	label N_ele = mesh.nCells(); 										//Read the number of Cells
	
	label Ns = matSnap.n(); 											//Number of snapshots
	
	//----------------------------------------------------------------------------------------------------//
	//Performing Singular Value Decomposition (SVD)
	DiagonalMatrix <scalar> singVal = SVD(matSnap,0).S(); 				//Singular values
	
	//----------------------------------------------------------------------------------------------------//
	//Sum of the Singular Values
	scalar sumSV = 0.0; 
	for (label i=0; i<Ns; i++)
	{
		sumSV += singVal[i];
	}
	
	//----------------------------------------------------------------------------------------------------//
	//Number of Dominant Modes [99.9999%]
	label Np = 0; 
	scalar var = 0.0;
	while (var <= 0.999999)
	{
		var += singVal[Np]/sumSV;
		Np += 1;
	}
	
	//----------------------------------------------------------------------------------------------------//
	//Formation of POD basis
	RectangularMatrix <scalar> vecLS = SVD(matSnap,0).U();
	RectangularMatrix <scalar> vecBasis(N_ele,Np);
	RectangularMatrix <scalar> vecBasisT(Np,N_ele);
	for (label i=0; i<N_ele; i++)
	{
		for (label j=0; j<Np; j++)
		{
			vecBasis[i][j] = vecLS[i][j];
		}
	}
	vecBasisT = vecBasis.T();
	
	//----------------------------------------------------------------------------------------------------//
	//SquareMatrix <scalar> matStiff(N_ele, 0.0); 						//Initialization of the Stiffness Matrix
	SquareMatrix <scalar> matStiff(N_ele, 0.0);
	RectangularMatrix <scalar> matTemp(N_ele,Np);
	RectangularMatrix <scalar> matProj(Np,Np);
	RectangularMatrix <scalar> matProjInv(Np,Np);
	
    //----------------------------------------------------------------------------------------------------//
    //Formation of the Stiffness Matrix
    while (runTime.loop())
    {
        fvScalarMatrix hEqn
		(
			fvm::laplacian(T, h)
			-
			fvm::ddt(b*Ss, h)
			-
			(b*Rc)
			+
			(b*Lk)
		);
        
        //Adding the Diagonal Coefficients
		for(label i=0; i<N_ele; i++) 
		{
			matStiff[i][i] = hEqn.diag()[i];
		}
		
		//Adding the Off-Diagonal Coefficients
		for(label faceID=0; faceID<hEqn.lduAddr().lowerAddr().size(); faceID++) 
		{
			label l = hEqn.lduAddr().lowerAddr()[faceID];
			label u = hEqn.lduAddr().upperAddr()[faceID];
			matStiff[l][u] = hEqn.upper()[faceID];
			matStiff[u][l] = hEqn.lower()[faceID];
		}
		
		//Adding contribution from Boundary Condition 
		forAll(h.boundaryField(), patchID) 
		{
			const fvPatch &patchType = h.boundaryField()[patchID].patch();
			forAll(patchType, faceID)
			{
				label cellID = patchType.faceCells()[faceID];
				matStiff[cellID][cellID] += hEqn.internalCoeffs()[patchID][faceID];
			}
		}
			
		//Stiffness Matrix X POD-Basis 
		matTemp = matStiff*vecBasis;
			
		//POD-Basis[Transpose] X matTemp
		matProj = vecBasisT*matTemp;
			
		//Inverse of Projected Stiffness Matrix
		matProjInv = SVDinv(matProj);
		
		break;
	}
	
	//----------------------------------------------------------------------------------------------------//
	std::clock_t endT= std::clock(); 									//End Time
	scalar simTime;
	simTime = (endT - startT)/ (double) CLOCKS_PER_SEC;
	
	//----------------------------------------------------------------------------------------------------//
	//Writing the Singular Values in a CSV file
	std::ofstream file1;
	file1.open("singularValues.csv");
	for (label i=0; i<Ns; i++)
	{
		file1 << singVal[i] << nl;
	}
	file1.close();
	
	//----------------------------------------------------------------------------------------------------//
	//Writing the POD basis in a CSV file
	std::ofstream file2;
	file2.open("PODBasisVectors.csv");
	for (label i=0; i<N_ele; i++)
	{
		for (label j=0; j<Np; j++)
		{
			file2 << vecBasis[i][j] << ",";
		}
		file2 << nl;
	}
	file2.close();
	
	//----------------------------------------------------------------------------------------------------//
	//Writing the Inverse of the Projected Stiffness Matrix in a CSV file
	std::ofstream file3;
	file3.open("InverseProjectedStiffnessMatrix.csv");
	for (label i=0; i<Np; i++)
	{
		for (label j=0; j<Np; j++)
		{
			file3 << matProjInv[i][j] << ",";
		}
		file3 << nl;
	}
	file3.close();
	
	//----------------------------------------------------------------------------------------------------//
	Info<< "\nCPU simulation time =" << simTime << nl << endl;
    Info<< "End\n" << endl;
    
	return 0;
}
