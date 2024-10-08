/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::radiation::collimatedbeamMixedFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a grey-diffuse condition for radiation
    intensity, \c I, for use with the finite-volume discrete-ordinates model
    (oqfvDOM), in which the radiation temperature is retrieved from the
    temperature field boundary condition.

    An external radiative heat flux can be added using \c qRadExt. If
    \c qRadExtDir is specified, this ray closest to this direction is used.
    Otherwise, the face normal is used as direction to set \c qRadExt.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            collimatedbeam;

        // Optional entries
        T               T;
        qRadExt         <scalar>;
        qRadExtDir      <vector>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                            | Type  | Reqd | Deflt
      type     | Type name: collimatedbeam      | word  | yes  | -
      T        | Name of temperature field              | word  | no   | T
      qRadExt  | Radiative external flux                | scalar | no  | Zero
      qRadExtDir | Radiative external flux direction    | vector | no  | Zero
    \endtable

    The inherited entries are elaborated in:
      - \link mixedFvPatchFields.H \endlink

See also
  - Foam::radiation::radiationModel
  - Foam::radiation::oqfvDOM
  - Foam::radiationCoupledBase
  - Foam::mixedFvPatchField

SourceFiles
    collimatedbeamMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_collimatedbeamMixedFvPatchScalarField_H
#define radiation_collimatedbeamMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
/*---------------------------------------------------------------------------*\
        Class collimatedbeamMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class collimatedbeamMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Name of temperature field
        word TName_;

        //- External radiative flux
        scalar qRadExt_;

        //- External radiative flux direction
        vector qRadExtDir_;
       
        //- irradiation                     //added
	scalar irradiation_;

	//- Beam Direction                      //added
	vector beamDirection_;  


public:

    //- Runtime type information
    TypeName("collimatedbeam");


    // Constructors

        //- Construct from patch and internal field
        collimatedbeamMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        collimatedbeamMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        //- collimatedbeamMixedFvPatchScalarField onto a new patch
        collimatedbeamMixedFvPatchScalarField
        (
            const collimatedbeamMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        collimatedbeamMixedFvPatchScalarField
        (
            const collimatedbeamMixedFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        collimatedbeamMixedFvPatchScalarField
        (
            const collimatedbeamMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member Functions

        // Access

            //- Return the temperature field name
            const word& TName() const noexcept
            {
                return TName_;
            }

            //- Return reference to the temperature field name to allow
            //- adjustment
            word& TName() noexcept
            {
                return TName_;
            }

    // const tmp<scalarField> temissivity();


        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
