/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 OpenCFD Ltd.
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
    Foam::radiation::lblBandAbsorptionEmission

Group
    grpRadiationAbsorptionEmissionSubModels

Description

    lblBandAbsorptionEmission radiation absorption/emission for solids.


SourceFiles
    lblBandAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef lblBandAbsorptionEmission_H
#define lblBandAbsorptionEmission_H

#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class lblBandAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class lblBandAbsorptionEmission
:
    public absorptionEmissionModel
{
public:

    // Public data

        //- Maximum number of bands
        static const int maxBands_ = 1000000;


private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //-  Absorption coefficients
        scalarList absCoeffs_;

        //-  Emissivity coefficients
        scalarList emiCoeffs_;

        //- Bands
        label nBands_;


public:

    //- Runtime type information
    TypeName("lblBandAbsorptionEmission");


    // Constructors

        //- Construct from components
        lblBandAbsorptionEmission
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~lblBandAbsorptionEmission();


    // Member Functions

        // Access
         //- Absorption coefficient (net)
                virtual tmp<volScalarField> a(const label bandI = 0) const;

            // Absorption coefficient

                //- Absorption coefficient
                tmp<volScalarField> aCont(const label bandI) const;


            // Emission coefficient

                //- Emission coefficient
                tmp<volScalarField> eCont(const label bandI) const;


            // Emission contribution

                //- Emission contribution
                tmp<volScalarField> ECont(const label bandI) const;


        inline bool isGrey() const
        {
            return false;
        }

        //- Number of bands
        inline label nBands() const
        {
            return nBands_;
        }
        //- Correct absorption coefficients
            virtual void correct
            (
                volScalarField& a,
                PtrList<volScalarField>& aj
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
