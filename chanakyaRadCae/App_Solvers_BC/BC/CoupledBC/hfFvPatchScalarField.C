/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "hfFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hfFvPatchScalarField::
hfFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    QrName_("qr"),
    thCond_(p.size(),0.0),    
    extFlux_(p.size(),0.0),
    unitVect_(pTraits<vector>::zero)    
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::hfFvPatchScalarField::
hfFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    QrName_(dict.lookup("qr")),
    thCond_("thCond", dict, p.size()),
    extFlux_("extFlux", dict, p.size()),
    unitVect_(dict.lookup("UnitVector"))     
{
    this->patchType() = dict.lookupOrDefault<word>("patchType", word::null);

    this->refValue() = 0.0;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(this->refValue());
    }

    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::hfFvPatchScalarField::
hfFvPatchScalarField
(
    const hfFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    QrName_(ptf.QrName_),
    thCond_(ptf.thCond_),
    extFlux_(ptf.thCond_),
    unitVect_(ptf.unitVect_)   
   
    
{}


Foam::hfFvPatchScalarField::
hfFvPatchScalarField
(
    const hfFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    QrName_(ptf.QrName_),
    thCond_(ptf.thCond_),
    extFlux_(ptf.thCond_),
    unitVect_(ptf.unitVect_) 
{}


Foam::hfFvPatchScalarField::
hfFvPatchScalarField
(
    const hfFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    QrName_(ptf.QrName_),
    thCond_(ptf.thCond_),
    extFlux_(ptf.thCond_),
    unitVect_(ptf.unitVect_) 
    
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hfFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalarField& Qr =patch().lookupPatchField<volScalarField, scalar>("qr");
   
    const vectorField n(patch().nf());
  // Info<<"the wall on which i am "<<patch().name()<<"The direction on the same wall is"<<patch().nf()<<endl;
   
     
      const scalarField unitFaceNormal(n & unitVect_);  
     
     // const scalarField nQr=Qr*unitFaceNormal;
     
   //Info<<"the direction changed vector normal"<<unitFaceNormal<<"The direction on the same wall is"<<nQr<<endl;
     
      valueFraction() =0.0;
      refValue() =0.0;
      refGrad()=(Qr+extFlux_)/thCond_;    
  
      
      mixedFvPatchScalarField::updateCoeffs();
}


void Foam::hfFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("qr") << QrName_ << token::END_STATEMENT << nl;  
    os.writeKeyword("thCond") <<thCond_ << token::END_STATEMENT << nl;  
    os.writeKeyword("extFlux") <<extFlux_ << token::END_STATEMENT << nl;  
    os.writeKeyword("UnitVector") <<unitVect_<< token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::hfFvPatchScalarField::operator=
(
    const fvPatchScalarField& ptf
)
{
    fvPatchScalarField::operator=
    (
        valueFraction()*refValue() + (1 - valueFraction())*ptf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        hfFvPatchScalarField
    );
}

// ************************************************************************* //
