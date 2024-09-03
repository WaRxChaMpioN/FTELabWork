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
    Foam::qinFvPatchScalarField

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
        type            ConvRad;
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
    qinFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef qinFvPatchScalarField_H
#define qinFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class qinFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class qinFvPatchScalarField
:
    public mixedFvPatchScalarField
{

protected:

    // Protected data

         //- free stream heat coefficient
        scalarField hinf_;

        //- free stream temperature
        scalarField Tinf_;

        //- conductivity of the medium
        scalarField Kappa_;
        
        //- Name of total radiation field from the boundary
        word QRName_;
        
         //- Name of incident radiation field from the boundary
        word QinName_;
        
        //- Name of emitted radiation field from the boundary
        word QemName_;
        
        //- define the unit Normal vector for patch
        vector unitVect_;


public:

    //- Runtime type information
    TypeName("ConvRadQin");


    // Constructors

        //- Construct from patch and internal field
        qinFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        qinFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  qinFvPatchScalarField onto a new patch
        qinFvPatchScalarField
        (
            const qinFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        qinFvPatchScalarField
        (
            const qinFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField >
            (
                new qinFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        qinFvPatchScalarField
        (
            const qinFvPatchScalarField&,
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
                new qinFvPatchScalarField(*this, iF)
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
            const word& QRName() const
            {
                return QRName_;
            }

            //- Return reference to the phaseFraction to allow adjustment
            word& QRName()
            {
                return QRName_;
            }
         
           //- Return the phaseFraction
            const word& QinName() const
            {
                return QinName_;
            }

            //- Return reference to the phaseFraction to allow adjustment
            word& QinName()
            {
                return QinName_;
            }
           
           
           //- Return the phaseFraction
            const word& QemName() const
            {
                return QemName_;
            }

            //- Return reference to the phaseFraction to allow adjustment
            word& QemName()
            {
                return QemName_;
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
