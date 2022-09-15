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
    Foam::readCSV[preComputation]

Group
    preComputation/confinedGWFlowROM/modFlowFOAM

Description
    Reads a CSV file line by line and stores in a Matrix.
\*-----------------------------------------------------------------------*/
#ifndef readCSV_H
#define readCSV_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "IFstream.H"
#include "OFstream.H"
#include "RectangularMatrix.H"
#include "List.H"

namespace Foam
{
	inline RectangularMatrix <scalar> readCSV(string filename)
	{
		std::vector<double> V;

		std::ifstream myfile(filename);
		string myline;
		string myentry;
		label nRow = 0;
		label nCol;
		while (getline(myfile, myline))
		{
			std::stringstream mystring(myline);
			nCol = 0;
			while (getline(mystring, myentry, ','))
			{
				V.push_back(std::stod(myentry));
				nCol++;
			}
			nRow++;
		}
	
		RectangularMatrix <scalar> M(nRow,nCol);
		for (label i=0; i<nRow; i++)
		{
			for (label j=0; j<nCol; j++)
			{
				M[i][j] = V[i*nCol+j];
			}
		}
		return M;
	}
}
#endif
