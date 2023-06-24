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
    Foam::pumpingConfinedBCFvPatchScalarField

Group
    pumpingConfinedBC/modFlowBoundaryConditions

Description
    This boundary condition calculates the normal gradient of pressure head 
    due to time constant pumping rate at the well-boundary for a confined 
    aquifer.
\*-----------------------------------------------------------------------*/
#include "pumpingConfinedBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

namespace Foam
{
	//Constructors
	
	pumpingConfinedBCFvPatchScalarField::pumpingConfinedBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(p, iF),
		qC_(0.0)
		{}
	/*-------------------------------------------------------------------------------------*/
	pumpingConfinedBCFvPatchScalarField::pumpingConfinedBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const dictionary& dict
	)
	:
		fixedGradientFvPatchScalarField(p, iF),
		qC_("qC", dict, p.size())
		{
			if (dict.found("gradient"))
			{
				gradient() = scalarField("gradient", dict, p.size());
				fixedGradientFvPatchScalarField::updateCoeffs();
				fixedGradientFvPatchScalarField::evaluate();
			}
			else
			{
				fvPatchField<scalar>::operator=(patchInternalField());
				gradient() = 0.0;
			}
		}
	/*-------------------------------------------------------------------------------------*/	
	pumpingConfinedBCFvPatchScalarField::pumpingConfinedBCFvPatchScalarField
	(
		const pumpingConfinedBCFvPatchScalarField& ptf,
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
		qC_(ptf.qC_, mapper)
		{}
	/*-------------------------------------------------------------------------------------*/
	pumpingConfinedBCFvPatchScalarField::pumpingConfinedBCFvPatchScalarField
	(
		const pumpingConfinedBCFvPatchScalarField& ptf
	)
	:
		fixedGradientFvPatchScalarField(ptf),
		qC_(ptf.qC_)
		{}
	/*-------------------------------------------------------------------------------------*/
	pumpingConfinedBCFvPatchScalarField::pumpingConfinedBCFvPatchScalarField
	(
		const pumpingConfinedBCFvPatchScalarField& ptf,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(ptf, iF),
		qC_(ptf.qC_)
		{}
	
	//Member Functions
	
	void pumpingConfinedBCFvPatchScalarField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}
		
		const fvPatchField<tensor>& T_B = patch().lookupPatchField<volTensorField, tensor>("T");
		
		const vectorField qC_B = qC_*(patch().nf());
                
        gradient() = -((qC_B & inv(T_B)) & patch().nf());
        
        fixedGradientFvPatchScalarField::updateCoeffs();
	}
	
	//Write
	void pumpingConfinedBCFvPatchScalarField::write(Ostream& os) const
	{
		fixedGradientFvPatchScalarField::write(os);
		qC_.writeEntry("qC", os);
		writeEntry("value", os);
	}
	
	makePatchTypeField(fvPatchScalarField, pumpingConfinedBCFvPatchScalarField);
			
} //End namespace Foam
