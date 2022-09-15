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
    Foam::reducedConfinedGWFlowFOAM

Group
    reducedConfinedGWFlowFOAM/confinedGWFlowROM/modFlowFOAM

Description
    Proper Orthogonal Decomposition (POD) based reduced-order modeling of 
    unsteady flow through a homogeneous/heterogeneous confined aquifer.
\*-----------------------------------------------------------------------*/
#include "fvCFD.H"
#include "fvOptions.H"
#include "IFstream.H"
#include "OFstream.H"
#include "readCSV.H"
#include "simpleMatrix.H"

int main(int argc, char *argv[])
{
	std::clock_t startT= std::clock(); 									//Start Time
	
    argList::addNote
    (
        "Reduced-order confined groundwater flow through a homogeneous/heterogeneous confined aquifer."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "RectangularMatrix.H"
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
	//Reading the POD-Basis from CSV files into a Matrix
	RectangularMatrix <scalar> vecBasis = readCSV("PODBasisVectors.csv");
	
	label N_ele = mesh.nCells(); 										//Read the number of Cells
	
	label Np = vecBasis.n(); 											//Number of POD-basis
	
	RectangularMatrix <scalar> vecBasisT(Np,N_ele);
	vecBasisT = vecBasis.T();
	
	//----------------------------------------------------------------------------------------------------//
	//Reading the Inverse of the Projected Stiffness Matrix from CSV files into a Matrix
	RectangularMatrix <scalar> matProjInv = readCSV("InverseProjectedStiffnessMatrix.csv");
		
	List<scalar> vecSource(N_ele); 										//Initializing the Source Vector
	List<scalar> vecProj(Np);
	List<scalar> coefTimeDep(Np); 										//Initializing the Time-Dependent Coefficient Vector
	List<scalar> hField(N_ele);
	
	//----------------------------------------------------------------------------------------------------//
	//Reduced-order model
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
		
		//Formation of the Source Vector
		for (label i=0; i<N_ele; i++)
		{
			vecSource[i] = hEqn.source()[i];
		}
		
		//Adding contribution from Boundary Condition
		forAll(h.boundaryField(), patchID) 
		{
			const fvPatch &patchType = h.boundaryField()[patchID].patch();
			forAll(patchType, faceID)
			{
				label cellID = patchType.faceCells()[faceID];
				vecSource[cellID] += hEqn.boundaryCoeffs()[patchID][faceID];
			}
		}
		
		//Projection of the Source Vector
		vecProj = vecBasisT*vecSource;
		
		//Solving for the Time-dependent Coefficients
		coefTimeDep = matProjInv*vecProj;
		
		//Reconstruction of the Hydraulic head Field
		hField = vecBasis*coefTimeDep;
		
		forAll(mesh.C(), cellID)
        {
			h[cellID] = hField[cellID];
		}
		
		h.correctBoundaryConditions();
		
		runTime.write();
	}
	
	//----------------------------------------------------------------------------------------------------//
	std::clock_t endT= std::clock(); 									//End Time
	scalar simTime;
	simTime = (endT - startT)/ (double) CLOCKS_PER_SEC;
	Info<< "\nCPU simulation time =" << simTime << nl << endl;
    Info<< "End\n" << endl;
    
	return 0;
}
