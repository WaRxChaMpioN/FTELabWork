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

#include "ConvRadFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ConvRadFvPatchScalarField::
ConvRadFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    hinf_(p.size(),0.0),
    Tinf_(p.size(),0.0),
    Kappa_(p.size(),0.0),
    QRName_("qr")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::ConvRadFvPatchScalarField::
ConvRadFvPatchScalarField
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
    QRName_(dict.lookup("qr"))
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


Foam::ConvRadFvPatchScalarField::
ConvRadFvPatchScalarField
(
    const ConvRadFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    hinf_(ptf.hinf_),
    Tinf_(ptf.Tinf_),
    Kappa_(ptf.Kappa_),
    QRName_(ptf.QRName_)
{}


Foam::ConvRadFvPatchScalarField::
ConvRadFvPatchScalarField
(
    const ConvRadFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf)   
{}


Foam::ConvRadFvPatchScalarField::
ConvRadFvPatchScalarField
(
    const ConvRadFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    hinf_(ptf.hinf_),
    Tinf_(ptf.Tinf_),
    Kappa_(ptf.Kappa_),
    QRName_(ptf.QRName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ConvRadFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // scalar rhor = 1000;
    // scalarField alphap1 = max(min(alphap, 1.0), 0.0);
    // valueFraction() = alphap1/(alphap1 + rhor*(1.0 - alphap1));
    
   /*  scalarField delta_=1/(patch().deltaCoeffs());
    
    scalarField Kbydelta_=Kappa_/delta_;
    
    scalarField beta_=Kappa_/(hinf_*delta_);
    
    this->valueFraction() = 1/(1+beta_);

    this->refValue() = Tinf_;
    
    this->refGrad() = 0.0;
  
   scalarField Tp =*this;
    
       // patch().lookupPatchField<volScalarField, scalar>(TName_);
    
    scalarField delta_=patch().deltaCoeffs();    
    scalarField Kbydelta_=Kappa_/delta_;
    
    scalarField ratio=hinf_/(hinf_+Kbydelta_);
    scalarField ratM=1-ratio;
    
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    
     volScalarField TW= mesh.lookupObject<volScalarField>("T"); //creating the object to hold the temperature field 
    //to hold the cell centre values for the patch which we are looping
     const labelList& Int_cnt=patch().faceCells();
     
     forAll(Int_cnt,faceI)
     {
     
      label celli=Int_cnt[faceI];
      
         
      Tp[faceI] = (TW[celli]);  
     
     }
    
    refValue()=(ratio*Tinf_)+(ratM*Tp);
    
    valueFraction()=1.0;
    */ 
    
    const scalarField& QR =
        patch().lookupPatchField<volScalarField, scalar>("qr");
    
    scalarField delta_=1/(patch().deltaCoeffs());
    
   // Info<<" The delta values are"<<delta_<<endl;
    
    scalarField Kbydelta_=Kappa_/delta_;
    
   // Info<<" The Kbydelta_ values are"<<Kbydelta_<<endl;
    
    scalarField beta_=(Kappa_/(hinf_*delta_));
    
    scalarField be_=(hinf_/(Kbydelta_+ hinf_));
    
    //Info<<" The be_ values are"<<be_<<endl;
    
    this->valueFraction() =(hinf_/(Kbydelta_+hinf_));

    this->refValue() = (Tinf_-(QR/hinf_));
   
  
    
    this->refGrad() = 0.0;
    
    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::ConvRadFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("hinf") << hinf_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tinf") << Tinf_ << token::END_STATEMENT << nl;
    os.writeKeyword("Kappa") << Kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("qr") << QRName_ << token::END_STATEMENT << nl; 
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::ConvRadFvPatchScalarField::operator=
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
        ConvRadFvPatchScalarField
    );
}

// ************************************************************************* //
