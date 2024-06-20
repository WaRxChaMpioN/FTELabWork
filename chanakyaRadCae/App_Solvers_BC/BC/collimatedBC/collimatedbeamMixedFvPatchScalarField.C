/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "collimatedbeamMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "boundaryRadiationProperties.H"

#include "oqfvDOM.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::collimatedbeamMixedFvPatchScalarField::
collimatedbeamMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    solarLoad_(false),
    irradiation_(0),
    beamDirection_(pTraits<vector>::zero)
    
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::radiation::collimatedbeamMixedFvPatchScalarField::
collimatedbeamMixedFvPatchScalarField
(
    const collimatedbeamMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    solarLoad_(ptf.solarLoad_),
    irradiation_(ptf.irradiation_),
    beamDirection_(ptf.beamDirection_)
{}


Foam::radiation::collimatedbeamMixedFvPatchScalarField::
collimatedbeamMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    solarLoad_(dict.lookupOrDefault<bool>("solarLoad", false)),
    irradiation_(readScalar(dict.lookup("Irradiation"))),
    beamDirection_(dict.lookup("BeamDirection"))
{
    if (dict.found("refValue"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        refValue() = 0.0;
        refGrad() = 0.0;
        valueFraction() = 1.0;

        fvPatchScalarField::operator=(refValue());
    }
}


Foam::radiation::collimatedbeamMixedFvPatchScalarField::
collimatedbeamMixedFvPatchScalarField
(
    const collimatedbeamMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    solarLoad_(ptf.solarLoad_),
    irradiation_(ptf.irradiation_),
    beamDirection_(ptf.beamDirection_)
{}


Foam::radiation::collimatedbeamMixedFvPatchScalarField::
collimatedbeamMixedFvPatchScalarField
(
    const collimatedbeamMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    solarLoad_(ptf.solarLoad_),
    irradiation_(ptf.irradiation_),
    beamDirection_(ptf.beamDirection_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::collimatedbeamMixedFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    const scalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    const radiationModel& radiation =
        db().lookupObject<radiationModel>("radiationProperties");

    const oqfvDOM& dom(refCast<const oqfvDOM>(radiation));

    label rayId = -1;
    label lambdaId = -1;
    dom.setRayIdLambdaId(internalField().name(), rayId, lambdaId);

    const label patchi = patch().index();

    if (dom.nLambda() != 1)
    {
        FatalErrorInFunction
            << " a grey boundary condition is used with a non-grey "
            << "absorption model" << nl << exit(FatalError);
    }

    scalarField& Iw = *this;

    const vectorField n(patch().nf());

    oqradiativeIntensityRay& ray =
        const_cast<oqradiativeIntensityRay&>(dom.IRay(rayId));

    const scalarField nAve(n & ray.dAve());

   
   ray.qr().boundaryFieldRef()[patchi] += Iw*nAve;
   
   scalar Omega=ray.omega();
   
    const boundaryRadiationProperties& boundaryRadiation =
        boundaryRadiationProperties::New(internalField().mesh());

    const tmp<scalarField> temissivity
    (
        boundaryRadiation.emissivity(patch().index())
    );

    const scalarField& emissivity = temissivity();

    scalarField& qem = ray.qem().boundaryFieldRef()[patchi];
    scalarField& qin = ray.qin().boundaryFieldRef()[patchi];

    const vector& myRayId = dom.IRay(rayId).d();

    // Use updated Ir while iterating over rays
    // avoids to used lagged qin
    scalarField Ir = dom.IRay(0).qin().boundaryField()[patchi];

    for (label rayI=1; rayI < dom.nRay(); rayI++)
    {
        Ir += dom.IRay(rayI).qin().boundaryField()[patchi];
     
    }
 
     
       
    if (solarLoad_)
    {
        Ir += patch().lookupPatchField<volScalarField,scalar>
        (
            radiation.externalRadHeatFieldName_
        );
    }


     forAll(Iw, faceI)
     {    
         
         if ((-n[faceI] & myRayId) > 0.0)
        {
                scalar diff1 = 0;
		scalar diff2 = 0;
		scalar diff3 = 0;
		if(myRayId[0] != 0.0)
		{
			diff1 = (beamDirection_[0] / myRayId[0]);
			if(diff1 > 0.99 && diff1 < 1.01)
			{
				diff1 = 1;
			}
		}
		else if(myRayId[0] == beamDirection_[0])
		{
			diff1 = 1;
		}
		if(myRayId[1] != 0.0)
		{
			//diff2 = round(beamDirection_[1] / myRayId[1]);
			diff2 = (beamDirection_[1] / myRayId[1]);
			if(diff2 > 0.99 && diff2 < 1.01)
			{
				diff2 = 1;
			}
		}
		else if(myRayId[1] == beamDirection_[1])
		{
			diff2 = 1;
		}
		if(myRayId[2] != 0.0)
		{
			//diff3 = round(beamDirection_[2] / myRayId[2]);
			diff3 = (beamDirection_[2] / myRayId[2]);
			if(diff3 > 0.99 && diff3 < 1.01)
			{
				diff3 = 1;
			}
		}
		else if(myRayId[2] == beamDirection_[2])
		{
			diff3 = 1;
		}
	
		if((diff1) == 1 && (diff2) == 1 && (diff3) == 1)
		{	    
		  refGrad()[faceI] = 0.0;
                   valueFraction()[faceI] = 1.0; 
		   refValue()[faceI] = fabs(irradiation_/Omega) +	   
		   (
                    Ir[faceI]*(scalar(1.0) - emissivity[faceI])
                  + emissivity[faceI]*physicoChemical::sigma.value()
                  * pow4(Tp[faceI])
                )/pi;
                 
	        qem[faceI] = refValue()[faceI]*nAve[faceI];
	    
		}
	
		else
		{
		 refGrad()[faceI] = 0.0;
                 valueFraction()[faceI] = 1.0;
		 refValue()[faceI] =
                (
                    Ir[faceI]*(scalar(1.0) - emissivity[faceI])
                  + emissivity[faceI]*physicoChemical::sigma.value()
                  * pow4(Tp[faceI])
                )/pi;
                
                qem[faceI] =refValue()[faceI]*nAve[faceI];
              
		}     
          
        }
      else
        {
            // direction into the wall
            valueFraction()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            refValue()[faceI] = 0.0; //not used

            // Incident heat flux on this ray direction
            qin[faceI] = (Iw[faceI]*nAve[faceI]);
        }
     
          
     }
   
  // Info<<"the qem values from the boundary"<<qem<<endl;
 // ray.qr().boundaryFieldRef()[patchi] = qin+qem;

 
 
    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::radiation::collimatedbeamMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("T", "T", TName_);
    os.writeKeyword("solarLoad") << solarLoad_ << token::END_STATEMENT << nl;
    os.writeKeyword("Irradiation") << irradiation_ << token::END_STATEMENT << nl;
    os.writeKeyword("BeamDirection") << beamDirection_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        collimatedbeamMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
