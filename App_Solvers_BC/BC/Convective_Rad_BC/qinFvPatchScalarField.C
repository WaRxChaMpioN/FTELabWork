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

#include "qinFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::qinFvPatchScalarField::
qinFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    hinf_(p.size(),0.0),
    Tinf_(p.size(),0.0),
    Kappa_(p.size(),0.0),
    QRName_("qr"),
    QinName_("qin"),
    QemName_("qem"),
    unitVect_(pTraits<vector>::zero)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::qinFvPatchScalarField::
qinFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    hinf_("hinf", dict, p.size()),
    Tinf_("Tinf", dict, p.size()),
    Kappa_("Kappa", dict, p.size()),
    QRName_(dict.lookup("qr")),
    QinName_(dict.lookup("qin")),
    QemName_(dict.lookup("qem")),
    unitVect_(dict.lookup("UnitVector"))
{
    this->patchType() = dict.lookupOrDefault<word>("patchType", word::null);

   // this->refValue() = pRefValue_;

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


Foam::qinFvPatchScalarField::
qinFvPatchScalarField
(
    const qinFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    hinf_(ptf.hinf_),
    Tinf_(ptf.Tinf_),
    Kappa_(ptf.Kappa_),
    QRName_(ptf.QRName_),
    QinName_(ptf.QinName_),
    QemName_(ptf.QemName_),
    unitVect_(ptf.unitVect_) 
{}


Foam::qinFvPatchScalarField::
qinFvPatchScalarField
(
    const qinFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf), 
    hinf_(ptf.hinf_),
    Tinf_(ptf.Tinf_),
    Kappa_(ptf.Kappa_),
    QRName_(ptf.QRName_),
    QinName_(ptf.QinName_),
    QemName_(ptf.QemName_),
    unitVect_(ptf.unitVect_)   
{}


Foam::qinFvPatchScalarField::
qinFvPatchScalarField
(
    const qinFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    hinf_(ptf.hinf_),
    Tinf_(ptf.Tinf_),
    Kappa_(ptf.Kappa_),
    QRName_(ptf.QRName_),
    QinName_(ptf.QinName_),
    QemName_(ptf.QemName_),
    unitVect_(ptf.unitVect_) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::qinFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

   
    const scalarField& QR =
        patch().lookupPatchField<volScalarField, scalar>("qr");
        
    //    const scalarField& Qin =
    //    patch().lookupPatchField<volScalarField, scalar>("qin");
        
    // const scalarField& Qem =
       // patch().lookupPatchField<volScalarField, scalar>("qem"); 
    
    scalarField delta_=1/(patch().deltaCoeffs());
    
   // Info<<" The delta values are"<<delta_<<endl;
    
    scalarField Kbydelta_=Kappa_/delta_;
    
    const vectorField n(patch().nf());
    
    const scalarField unitFaceNormal(n & unitVect_);  
    
    const scalarField nQr=QR*unitFaceNormal;
    
    scalarField qrbyhinf_=QR/hinf_;
    
    
    
    this->valueFraction() =(hinf_/(Kbydelta_+hinf_));

    this->refValue() = (Tinf_+qrbyhinf_);
    
    this->refGrad() = 0.0;
    
    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::qinFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("hinf") << hinf_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tinf") << Tinf_ << token::END_STATEMENT << nl;
    os.writeKeyword("Kappa") << Kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("qr") << QRName_ << token::END_STATEMENT << nl; 
    os.writeKeyword("qin") << QinName_ << token::END_STATEMENT << nl; 
    os.writeKeyword("qem") << QemName_ << token::END_STATEMENT << nl; 
    os.writeKeyword("UnitVector") <<unitVect_<< token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::qinFvPatchScalarField::operator=
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
        qinFvPatchScalarField
    );
}

// ************************************************************************* //
