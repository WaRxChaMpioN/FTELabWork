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

#include "collimatedFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::collimatedFluxFvPatchScalarField::
collimatedFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    QrName_("qr"),
    QemName_("qem"),
    QinName_("qin"),
    thCond_(p.size(),0.0),    
    extFlux_(p.size(),0.0)    
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::collimatedFluxFvPatchScalarField::
collimatedFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    QrName_(dict.lookup("qr")),
    QemName_(dict.lookup("qem")),
    QinName_(dict.lookup("qin")),
    thCond_("thCond", dict, p.size()),
    extFlux_("extFlux", dict, p.size())
    
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


Foam::collimatedFluxFvPatchScalarField::
collimatedFluxFvPatchScalarField
(
    const collimatedFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    QrName_(ptf.QrName_),
    QemName_(ptf.QemName_),
    QinName_(ptf.QinName_),
    thCond_(ptf.thCond_),
    extFlux_(ptf.thCond_)
    
{}


Foam::collimatedFluxFvPatchScalarField::
collimatedFluxFvPatchScalarField
(
    const collimatedFluxFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    QrName_(ptf.QrName_),
    QemName_(ptf.QemName_),
    QinName_(ptf.QinName_),
    thCond_(ptf.thCond_),
    extFlux_(ptf.thCond_)
{}


Foam::collimatedFluxFvPatchScalarField::
collimatedFluxFvPatchScalarField
(
    const collimatedFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    QrName_(ptf.QrName_),
    QemName_(ptf.QemName_),
    QinName_(ptf.QinName_),
    thCond_(ptf.thCond_),
    extFlux_(ptf.thCond_)
    
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::collimatedFluxFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalarField& Qr =
        patch().lookupPatchField<volScalarField, scalar>("qr");
        
    const scalarField& Qem =  patch().lookupPatchField<volScalarField, scalar>("qem");
    
    const scalarField& Qin =  patch().lookupPatchField<volScalarField, scalar>("qin");

   //  scalarField& Tw = *this;
     
      valueFraction() =0.0;
      refValue() =0.0;
      refGrad()=(extFlux_-Qr)/thCond_;    
  
      
      mixedFvPatchScalarField::updateCoeffs();
}


void Foam::collimatedFluxFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("qr") << QrName_ << token::END_STATEMENT << nl;  
    os.writeKeyword("qem") << QemName_ << token::END_STATEMENT << nl; 
    os.writeKeyword("qin") << QinName_ << token::END_STATEMENT << nl; 
    os.writeKeyword("thCond") <<thCond_ << token::END_STATEMENT << nl;  
    os.writeKeyword("extFlux") <<extFlux_ << token::END_STATEMENT << nl;  
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::collimatedFluxFvPatchScalarField::operator=
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
        collimatedFluxFvPatchScalarField
    );
}

// ************************************************************************* //
