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
/*---------------------------------------------------------------------------*\
Class
    Foam::wellConfinedBCFvPatchScalarField

Group
    wellConfinedBC/modFlowBoundaryConditions

Description
    This boundary condition calculates the drawdown at the well and surrounding 
    locations for unsteady groundwater flow through a confined aquifer.
\*-----------------------------------------------------------------------*/
#include "wellConfinedBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "besselK0.H"
#include "besselK1.H"
#include "stehfestCoefficients.H"

namespace Foam
{
	//Constructors
	
	wellConfinedBCFvPatchScalarField::wellConfinedBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedValueFvPatchScalarField(p, iF),
		Kr_(0.0),
		Kz_(0.0),
		Ss_(0.0),
		Qp_(0.0),
		rW_(0.0),
		r_(0.0),
		B_(0.0),
		d_(0.0),
		l_(0.0),
		dS_(0.0),
		Ksw_(0.0),
		rWC_(0.0)
		{}
	/*--------------------------------------------------------------------------*/
	wellConfinedBCFvPatchScalarField::wellConfinedBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const dictionary& dict
	)
	:
		fixedValueFvPatchScalarField(p, iF),
		Kr_(readScalar(dict.lookup("Kr"))),
		Kz_(readScalar(dict.lookup("Kz"))),
		Ss_(readScalar(dict.lookup("Ss"))),
		Qp_(readScalar(dict.lookup("Qp"))),
		rW_(readScalar(dict.lookup("rW"))),
		r_(readScalar(dict.lookup("r"))),
		B_(readScalar(dict.lookup("B"))),
		d_(readScalar(dict.lookup("d"))),
		l_(readScalar(dict.lookup("l"))),
		dS_(readScalar(dict.lookup("dS"))),
		Ksw_(readScalar(dict.lookup("Ksw"))),
		rWC_(readScalar(dict.lookup("rWC")))
		{
			if (dict.found("value"))
			{
				fvPatchField<scalar>::operator = (scalarField("value", dict, p.size()));
				fixedValueFvPatchScalarField::updateCoeffs();
				fixedValueFvPatchScalarField::evaluate();
			}
			else
			{
				fvPatchField<scalar>::operator=(patchInternalField());
			}
		}
	/*--------------------------------------------------------------------------*/
	wellConfinedBCFvPatchScalarField::wellConfinedBCFvPatchScalarField
	(
		const wellConfinedBCFvPatchScalarField& ptf,
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedValueFvPatchScalarField(ptf, p, iF, mapper),
		Kr_(ptf.Kr_),
		Kz_(ptf.Kz_),
		Ss_(ptf.Ss_),
		Qp_(ptf.Qp_),
		rW_(ptf.rW_),
		r_(ptf.r_),
		B_(ptf.B_),
		d_(ptf.d_),
		l_(ptf.l_),
		dS_(ptf.dS_),
		Ksw_(ptf.Ksw_),
		rWC_(ptf.rWC_)
		{}
	/*--------------------------------------------------------------------------*/
	wellConfinedBCFvPatchScalarField::wellConfinedBCFvPatchScalarField
	(
		const wellConfinedBCFvPatchScalarField& wucpsf
	)
	:
		fixedValueFvPatchScalarField(wucpsf),
		Kr_(wucpsf.Kr_),
		Kz_(wucpsf.Kz_),
		Ss_(wucpsf.Ss_),
		Qp_(wucpsf.Qp_),
		rW_(wucpsf.rW_),
		r_(wucpsf.r_),
		B_(wucpsf.B_),
		d_(wucpsf.d_),
		l_(wucpsf.l_),
		dS_(wucpsf.dS_),
		Ksw_(wucpsf.Ksw_),
		rWC_(wucpsf.rWC_)
		{}
	/*--------------------------------------------------------------------------*/
	wellConfinedBCFvPatchScalarField::wellConfinedBCFvPatchScalarField
	(
		const wellConfinedBCFvPatchScalarField& wucpsf,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedValueFvPatchScalarField(wucpsf, iF),
		Kr_(wucpsf.Kr_),
		Kz_(wucpsf.Kz_),
		Ss_(wucpsf.Ss_),
		Qp_(wucpsf.Qp_),
		rW_(wucpsf.rW_),
		r_(wucpsf.r_),
		B_(wucpsf.B_),
		d_(wucpsf.d_),
		l_(wucpsf.l_),
		dS_(wucpsf.dS_),
		Ksw_(wucpsf.Ksw_),
		rWC_(wucpsf.rWC_)
		{}
	/*--------------------------------------------------------------------------*/
	//Member Functions
	/*--------------------------------------------------------------------------*/
	void wellConfinedBCFvPatchScalarField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}
		
		constexpr scalar Pi = 3.141592654;
		constexpr scalar maxQ = 708.0;
		constexpr label Ns = 8;
		constexpr label Nc_min = 4;
		constexpr label Nc_max = 30;
		constexpr scalar zD = 0.5;
				
		const fvPatchField<scalar>& h = patch().lookupPatchField<volScalarField, scalar>("h");
		
		scalar t = this->db().time().value();
		
		//Dimensionless variables
		scalar rD = r_/rW_;									
		scalar rWD = rW_/B_;					  				
		scalar dD = d_/B_;					  				
		scalar lD = l_/B_;									
		
		scalar betaCD = (Kz_/Kr_)*sqr(rWD);
		scalar betaC = betaCD*sqr(rD);
		scalar SW = (Kr_/Ksw_)*(dS_/rW_);						//Dimensionless well-bore skin parameter
		scalar WD = sqr(rWC_)/(2.0*Ss_*(l_ - d_)*sqr(rW_)); 	//Dimensionless well-bore storage parameter 
		scalar tD = Kr_*t/(Ss_*sqr(rW_));						//Dimensionless time
		
		//Number of terms for the finite summation of drawdown
		label Nc = Nc_max*(pow(2.0,(-log10(betaC)-2.0)));		
		if (Nc < Nc_min)
		{
			Nc = Nc_min;
		}
		else if (Nc > Nc_max)
		{
			Nc = Nc_max;
		}
		
		scalar saWD = 0.0;
		scalarField saOD = 0.0*h;
		label ct = 1;
		
		for (label i = 1; i <= Ns; i++)
		{
			scalar p = ct*log(2.0)/tD;
			
			//Calculation of the Laplacian of the drwadowns at well and observation locations
			scalar A0 = 0.0;
			scalar A = 0.0;
			scalarField E0 = 0.0*h;
			scalarField E = 0.0*h;
			scalar eps = Pi;
			
			scalar q01 = sqrt(p);
			if (q01 > maxQ)
			{
				q01 = maxQ;
			}
			scalar q02 = q01*rD;
			if (q02 > maxQ)
			{
				q02 = maxQ;
			}
			
			A0 = besselK0(q01)*(lD - dD)/(q01*besselK1(q01));
			E0 = besselK0(q02)*(lD - dD)/(q01*besselK1(q01));
						
			for(label j = 1; j <= Nc; j++)
			{
				scalar qN1 = sqrt(p + sqr(eps)*betaCD);
				if (qN1 > maxQ)
				{
					qN1 = maxQ;
				}
				scalar qN2 = qN1*rD;
				if (qN2 > maxQ)
				{
					qN2 = maxQ;
				}
								
				scalar coef1 = sin(eps*(1.0 - dD)) - sin(eps*(1.0 - lD));
				scalar coef2 = eps*qN1*besselK1(qN1);
				
				A = A + (2.0/(lD - dD))*(besselK0(qN1)*sqr(coef1)/(eps*coef2));
				E = E + 2.0*(besselK0(qN2)*cos(eps*zD)*coef1/coef2);
				
				eps = eps + Pi; 
			} 
						
			scalar lap_saWD = (A0 + A + SW)/(p*(1.0 + WD*p*(A0 + A + SW)));
			scalarField lap_saOD = (E0 + E)/(p*(1.0 + WD*p*(A0 + A + SW)));
			
			scalar wt = stehfestCoef(ct, Ns);

			saWD = saWD + wt*lap_saWD;
			saOD = saOD + wt*lap_saOD;
			
			ct = ct + 1;
		} 
		
		saWD = (log(2.0)/tD)*(2.0/(lD - dD))*saWD;
		saOD = (log(2.0)/tD)*(2.0/(lD - dD))*saOD;
	
		scalar hW = B_ - saWD*Qp_/(4.0*Pi*Kr_*B_);
		scalarField hO = B_ - saOD*Qp_/(4.0*Pi*Kr_*B_);
		
		Info << "\nHydraulic head at the well =" << hW << endl;
		
		operator == (hO);

		fixedValueFvPatchScalarField::updateCoeffs();	
		
	}
	
	//Write
	void wellConfinedBCFvPatchScalarField::write(Ostream& os) const
	{
		fvPatchScalarField::write(os);
		os.writeEntry("Kr", Kr_);
		os.writeEntry("Kz", Kz_);
		os.writeEntry("Ss", Ss_);
		os.writeEntry("Qp", Qp_);
		os.writeEntry("rW", rW_);
		os.writeEntry("r", r_);
		os.writeEntry("B", B_);
		os.writeEntry("d", d_);
		os.writeEntry("l", l_);
		os.writeEntry("dS", dS_);
		os.writeEntry("Ksw", Ksw_);
		os.writeEntry("rWC", rWC_);
		writeEntry("value", os);
	}
	
	makePatchTypeField(fvPatchScalarField, wellConfinedBCFvPatchScalarField);
			
} //End namespace Foam
