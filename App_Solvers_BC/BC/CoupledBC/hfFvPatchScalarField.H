/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::hfFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a phase-based hydrostatic pressure
    condition, calculated as:

        \f[
            p_{hyd} = p_{ref} + \rho g (x - x_{ref})
        \f]

    where
    \vartable
        p_{hyd} | hyrostatic pressure [Pa]
        p_{ref} | reference pressure [Pa]
        x_{ref} | reference point in Cartesian co-ordinates
        \rho    | density (assumed uniform)
        g       | acceleration due to gravity [m/s2]
    \endtable

    The values are assigned according to the phase-fraction field:
    - 1: apply \$fp_{hyd}\$f
    - 0: apply a zero-gradient condition

Usage
    \table
        Property      | Description                 | Required | Default value
        phaseFraction | phase-fraction field name   | no       | alpha
        rho           | density field name          | no       | rho
        pRefValue     | reference pressure [Pa]     | yes      |
        pRefPoint     | reference pressure location | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            hf;
        phaseFraction   alpha1;
        rho             rho;
        pRefValue       1e5;
        pRefPoint       (0 0 0);
        value           uniform 0; // optional initial value
    }
    \endverbatim

See also
    Foam::mixedFvPatchScalarField

SourceFiles
    hfFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef hfFvPatchScalarField_H
#define hfFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class hfFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class hfFvPatchScalarField
:
    public mixedFvPatchScalarField
{

protected:

    // Protected data

        //- Name of phase-fraction field
        word QrName_;
        
   //- conductivity of the medium
        scalarField thCond_;
        
   //- external flux applied at the wall
        scalarField extFlux_;
        
    //- define the unit Normal vector for patch
        vector unitVect_;


public:

    //- Runtime type information
    TypeName("hf");


    // Constructors

        //- Construct from patch and internal field
        hfFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        hfFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  hfFvPatchScalarField onto a new patch
        hfFvPatchScalarField
        (
            const hfFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        hfFvPatchScalarField
        (
            const hfFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField >
            (
                new hfFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        hfFvPatchScalarField
        (
            const hfFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new hfFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Attributes

            //- Return true: this patch field is altered by assignment
            virtual bool assignable() const
            {
                return true;
            }


        // Access

            //- Return the phaseFraction
            const word& QrName() const
            {
                return QrName_;
            }

            //- Return reference to the phaseFraction to allow adjustment
            word& QrName()
            {
                return QrName_;
            }

          

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchScalarField& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
