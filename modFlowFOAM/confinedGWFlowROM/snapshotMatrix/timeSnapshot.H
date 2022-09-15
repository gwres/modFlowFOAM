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
    Foam::timeSnapshot[snapshotMatrix]

Group
    snapshotMatrix/confinedGWFlowROM/modFlowFOAM

Description
    Calculates the snapshot time-set based on an exponential function.
\*-----------------------------------------------------------------------*/
#ifndef timeSnapshot_H
#define timeSnapshot_H

#include "List.H"

namespace Foam
{
	inline List<scalar> dtSnap(const label Ns, const scalar c, const scalar Ts)
	{
		List<scalar> timeSet(Ns);
		List<scalar> dt(Ns-1);
		scalar A = 1.0/(Ns-1);
		scalar U = 0; 
		for (label i=0; i<Ns; i++)
		{
			timeSet[i] = Ts*(1.0 - tanh(c*(1 - U)));
			U += A;
		}
		timeSet[Ns-1] = Ts;
		for (label j=0; j<(Ns-1); j++)
		{
			dt[j] = timeSet[j+1] - timeSet[j];
		}
		return dt;
	}
} //End Namespace Foam

#endif
