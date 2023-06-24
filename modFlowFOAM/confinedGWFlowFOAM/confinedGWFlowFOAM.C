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
    Foam::confinedGWFlowFOAM

Group
    confinedGWFlowFOAM/modFlowFOAM

Description
    Modeling of unsteady state groundwater flow through a homogeneous/
    heterogeneous confined aquifer
\*-----------------------------------------------------------------------*/
#include "fvCFD.H"

int main(int argc, char *argv[])
{
    std::clock_t startT= std::clock(); 									//Start Time
    
    argList::addNote
    (
        "Confined groundwater flow through a homogeneous/heterogeneous confined aquifer."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    
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
	//Time-Loop
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

        runTime.write();
        
        Info << "-------------------------------------------------------" << endl;
    }

	//----------------------------------------------------------------------------------------------------//
	std::clock_t endT= std::clock(); 									//End Time
	
	scalar simTime;
	simTime = (endT - startT)/ (double) CLOCKS_PER_SEC;
	Info<< "\nCPU simulation time =" << simTime << nl << endl;
    
    Info<< "End\n" << endl;

    return 0;
}
