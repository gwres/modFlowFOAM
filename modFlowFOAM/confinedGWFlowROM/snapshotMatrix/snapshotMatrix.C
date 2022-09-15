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
    Foam::snapshotMatrix

Group
    snapshotMatrix/confinedGWFlowROM/modFlowFOAM

Description
    Formation of snapshot matrix for reduced-order modeling of unsteady 
    flow through a homogeneous/heterogeneous confined aquifer.
\*-----------------------------------------------------------------------*/
#include "fvCFD.H"
#include "fvOptions.H"
#include "IFstream.H"
#include "OFstream.H"
#include "timeSnapshot.H"
#include "SVD.H"

int main(int argc, char *argv[])
{
    std::clock_t startT= std::clock(); 									//Start Time
    
    argList::addNote
    (
        "Formation of Proper Orthogonal Decomposition basis"
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
	//Calculation of variable DeltaT for Snapshot Formation
    List<scalar> dt(Ns-1);
    dt = dtSnap(Ns, c, Ts);
    
    //----------------------------------------------------------------------------------------------------//
    //Calculating the Snapshots for unsteady state groundwater flow
    label N_ele = mesh.nCells(); 										//Read the number of Cells
     
    RectangularMatrix <scalar> matSnap(N_ele,Ns); 						//Declare the Snapshot Matrix
    
    forAll(h.internalField(), cellID)
    {
		matSnap[cellID][0] = h.internalField()[cellID];
	}
 
	//----------------------------------------------------------------------------------------------------//
	//Time-Loop
	label count = 1;
	runTime.setDeltaT(dt[0]);
    while (runTime.loop())
    {
        Info<< "\nTime = " << runTime.timeName() << nl << endl;
        
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
            
        hEqn.solve();

        forAll(h.internalField(), cellID)
        {
			matSnap[cellID][count] = h.internalField()[cellID];
		}
		
		runTime.write();
		
		runTime.setDeltaT(dt[count]);
		
		count += 1;
		
		Info << "-------------------------------------------------------" << endl;
		
		if (count == Ns)
		{
			break;
		}   
    }
	
	//----------------------------------------------------------------------------------------------------//
	std::clock_t endT = std::clock(); 									//End Time
	scalar simTime;
	simTime = (endT - startT)/ (double) CLOCKS_PER_SEC;
	
	//----------------------------------------------------------------------------------------------------//
	//Writing the snapshot matrix in a CSV file
	std::ofstream file1;
	file1.open("snapshotMatrix.csv");
	for (label i=0; i<N_ele; i++)
	{
		for (label j=0; j<Ns; j++)
		{
			file1 << matSnap[i][j] << ",";
		}
		file1 << nl;
	}
	file1.close();
	
	//----------------------------------------------------------------------------------------------------//
	//Writing the Snapshot Time-set in a CSV file
	std::ofstream file2;
	file2.open("snapshotTimeSet.csv");
	scalar start_time = 0;
	for (label i=0; i<Ns-1; i++)
	{
		start_time += dt[i];
		file2 << start_time << nl;
	}
	file2.close();
	
	//----------------------------------------------------------------------------------------------------//
	Info<< "\nCPU simulation time =" << simTime << nl << endl;
    Info<< "End\n" << endl;

    return 0;
}
