/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
// #include "boundaryRadiationProperties.C"

#include "oqfvDOM.H"
#include "constants.H"
#include "unitConversion.H"

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
    qRadExt_(0),
    qRadExtDir_(Zero),
    irradiation_(0),                        //added
    beamDirection_(pTraits<vector>::zero)  //added
{
    refValue() = Zero;
    refGrad() = Zero;
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
    qRadExt_(ptf.qRadExt_),
    qRadExtDir_(ptf.qRadExtDir_),
    irradiation_(ptf.irradiation_),             //added
    beamDirection_(ptf.beamDirection_)          //added
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
    TName_(dict.getOrDefault<word>("T", "T")),
    qRadExt_(dict.getOrDefault<scalar>("qRadExt", 0)),
    qRadExtDir_(dict.getOrDefault<vector>("qRadExtDir", Zero)),
    irradiation_(readScalar(dict.lookup("Irradiation"))),               //added
    beamDirection_(dict.lookup("BeamDirection"))                        //added
{
    if (this->readMixedEntries(dict))
    {
        // Full restart
        this->readValueEntry(dict, IOobjectOption::MUST_READ);
    }
    else
    {
        refValue() = Zero;
        refGrad() = Zero;
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
    qRadExt_(ptf.qRadExt_),
    qRadExtDir_(ptf.qRadExtDir_),
    irradiation_(ptf.irradiation_),                 //added
    beamDirection_(ptf.beamDirection_)              //added
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
    qRadExt_(ptf.qRadExt_),
    qRadExtDir_(ptf.qRadExtDir_),
    irradiation_(ptf.irradiation_),                 //added
    beamDirection_(ptf.beamDirection_)              //added
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Foam::tmp<Foam::scalarField> 
// Foam::radiation::collimatedbeamMixedFvPatchScalarField::temissivity()
// {
//     return boundaryRadiation.emissivity(patch().index(), 0, nullptr, &Tp)
// }
Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::emissivity
(
    const label patchi,
    const label bandi,
    const vectorField* incomingDirection,
    const scalarField* T
) const
{
    if (radBoundaryPropertiesPtrList_.set(patchi))
    {
        return radBoundaryPropertiesPtrList_[patchi].e
        (
            bandi,
            incomingDirection,
            T
        );
    }

    FatalErrorInFunction
        << "Patch : " << mesh().boundaryMesh()[patchi].name()
        << " is not found in the boundaryRadiationProperties. "
        << "Please add it"
        << exit(FatalError);

    return nullptr;
}

void Foam::radiation::collimatedbeamMixedFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    const auto& Tp = patch().lookupPatchField<volScalarField>(TName_);

    const oqfvDOM& dom = db().lookupObject<oqfvDOM>("radiationProperties");

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

   scalar Omega=ray.omega();                                //added

    const boundaryRadiationProperties& boundaryRadiation =
        boundaryRadiationProperties::New(internalField().mesh());

    const tmp<scalarField> temissivity
    (
        boundaryRadiation.emissivity(patch().index(), 0, nullptr, &Tp)
    );

    const scalarField& emissivity = temissivity();
    // const scalarField& emissivity = boundaryRadiation.emissivity(patch().index(), 0, nullptr, &Tp);

    // const tmp<scalarField> ttransmissivity
    // (
    //     boundaryRadiation.transmissivity(patch().index(), 0, nullptr, &Tp)
    // );

    // const scalarField& transmissivity = ttransmissivity();

    scalarField& qem = ray.qem().boundaryFieldRef()[patchi];
    scalarField& qin = ray.qin().boundaryFieldRef()[patchi];

    const vector& myRayId = dom.IRay(rayId).d();

    scalarField Ir(patch().size(), Zero);
    forAll(Iw, facei)
    {
        for (label rayi=0; rayi < dom.nRay(); rayi++)
        {
            const vector& d = dom.IRay(rayi).d();

            if ((-n[facei] & d) < 0.0)
            {
                // q into the wall
                const scalarField& IFace =
                    dom.IRay(rayi).ILambda(lambdaId).boundaryField()[patchi];

                const vector& rayDave = dom.IRay(rayi).dAve();
                Ir[facei] += IFace[facei]*(n[facei] & rayDave);
            }
        }
    }

    if (dom.useSolarLoad())
    {
        // Looking for primary heat flux single band
        Ir += patch().lookupPatchField<volScalarField>
        (
            dom.primaryFluxName_ + "_0"
        );

        word qSecName = dom.relfectedFluxName_ + "_0";

        if (this->db().foundObject<volScalarField>(qSecName))
        {
             const volScalarField& qSec =
                this->db().lookupObject<volScalarField>(qSecName);

            Ir += qSec.boundaryField()[patch().index()];
        }
    }

    scalarField Iexternal(this->size(), 0.0);

    if (dom.useExternalBeam())
    {
        const vector sunDir = dom.solarCalc().direction();
        const scalar directSolarRad = dom.solarCalc().directSolarRad();

        //label nRaysBeam = dom.nRaysBeam();
        label SunRayId(-1);
        scalar maxSunRay = -GREAT;

        // Looking for the ray closest to the Sun direction
        for (label rayI=0; rayI < dom.nRay(); rayI++)
        {
            const vector& iD = dom.IRay(rayI).d();
            scalar dir = sunDir & iD;
            if (dir > maxSunRay)
            {
                maxSunRay = dir;
                SunRayId = rayI;
            }
        }

        if (rayId == SunRayId)
        {
            const scalarField nAve(n & dom.IRay(rayId).dAve());
            forAll(Iexternal, faceI)
            {
                Iexternal[faceI] = directSolarRad/mag(dom.IRay(rayId).dAve());
            }
        }
    }

    scalarField Isource(this->size(), 0.0);

    if (qRadExt_ > 0)
    {
        if (mag(qRadExtDir_) > 0)
        {
            label rayqoId = -1;
            scalar maxRay = -GREAT;

            // Looking for the ray closest to the Sun direction
            for (label rayI = 0; rayI < dom.nRay(); ++rayI)
            {
                const vector& iD = dom.IRay(rayI).d();
                const scalar dir = qRadExtDir_ & iD;

                if (dir > maxRay)
                {
                    maxRay = dir;
                    rayqoId = rayI;
                }
            }

            if (rayId == rayqoId)
            {
                forAll(Isource, faceI)
                {
                    Isource[faceI] += qRadExt_/mag(dom.IRay(rayId).dAve());
                }
            }
        }
        else
        {
            forAll(Iw, faceI)
            {
                label rayqoId = -1;
                scalar maxRay = -GREAT;

                // Looking for the ray closest to the Sun direction
                for (label rayI = 0; rayI < dom.nRay(); ++rayI)
                {
                    const vector& iD = dom.IRay(rayI).d();
                    const scalar dir = -n[faceI] & iD;

                    if (dir > maxRay)
                    {
                        maxRay = dir;
                        rayqoId = rayI;
                    }
                }

                if (rayId == rayqoId)
                {
                    Isource[faceI] += qRadExt_/mag(dom.IRay(rayId).dAve());
                }
            }
        }
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
            qin[faceI] = Iw[faceI]*nAve[faceI];
        }
    }

    UPstream::msgType(oldTag);  // Restore tag

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::radiation::collimatedbeamMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("T", "T", TName_);
    os.writeEntryIfDifferent<scalar>("qRadExt", Zero, qRadExt_);
    os.writeEntryIfDifferent<vector>("qRadExtDir", Zero, qRadExtDir_);
    os.writeKeyword("BeamDirection") << beamDirection_ << token::END_STATEMENT << nl;
    os.writeKeyword("Irradiation") << irradiation_ << token::END_STATEMENT << nl;
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
