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
    Foam::createFields[confinedGWFlowFOAM]

Group
    confinedGWFlowFOAM/modFlowFOAM

Description
    Unsteady state groundwater flow through a homogeneous/heterogeneous 
    confined aquifer.
\*-----------------------------------------------------------------------*/
Info<< "Reading HYDRAULIC HEAD field h\n" << endl;

volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*-------------------------------------------------------------*/
Info<< "Reading SATURATED HYDRAULIC CONDUCTIVITY field Ks\n" << endl;

volTensorField Ks
(
    IOobject
    (
        "Ks",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*-------------------------------------------------------------*/
Info<< "Reading TRANSMISSIVITY field T\n" << endl;

volTensorField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*-------------------------------------------------------------*/
Info<< "Reading SPECIFIC STORAGE field Ss\n" << endl;

volScalarField Ss
(
    IOobject
    (
        "Ss",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*-------------------------------------------------------------*/
Info<< "Reading DEPTH field b\n" << endl;

volScalarField b
(
    IOobject
    (
        "b",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*-------------------------------------------------------------*/
Info<< "Reading RECHARGE TERM field Rc\n" << endl;

volScalarField Rc
(
    IOobject
    (
        "Rc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*-------------------------------------------------------------*/
Info<< "Reading LEAKAGE TERM field Lk\n" << endl;

volScalarField Lk
(
    IOobject
    (
        "Lk",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*-------------------------------------------------------------*/
// Reads user-defined data from modelDict in system directory

IOdictionary modelDict
(
	IOobject
    (
		"modelDict",
		runTime.system(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
    )
);
/*------------------------------------------------------------------------------*/
