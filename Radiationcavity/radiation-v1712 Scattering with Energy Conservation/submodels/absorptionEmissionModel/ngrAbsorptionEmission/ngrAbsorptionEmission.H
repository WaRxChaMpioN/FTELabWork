/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::radiation::ngrAbsorptionEmission

Group
    grpRadiationAbsorptionEmissionSubModels

Description

    ngrAbsorptionEmission radiation absorption/emission for solids.


SourceFiles
    ngrAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef ngrAbsorptionEmission_H
#define ngrAbsorptionEmission_H
#include "absorptionCoeffs.H"
#include "absorptionEmissionModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class ngrAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class ngrAbsorptionEmission
:
    public absorptionEmissionModel
{
public:

    // Public data

        //- Maximum number of bands
        static const int maxBands_ =32;
		
				


private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //-  Absorption coefficients
        scalarList absCoeffs_;

        //-  Emissivity coefficients
        scalarList emiCoeffs_;

        //- Number of weights
        scalarList weights_;

        //- Bands
        label nBands_;

        //- Total number of weights used
        label nWeights_;

        //- Bands
        FixedList<Vector2D<scalar>, maxBands_> iBands_;

        // Total wave length covered by the bands
        scalar totalWaveLength_;


        //Total set of moleFractions datas available in dataBase
        static const int nMoleFractions=12;
			
        // Minimum Temperature at which database is available
		static const int minTemp=300;

        // Maximum Temperature at which database is available
		static const int maxTemp=2100;

        // Temperature step for which database is available
        static const int stepTemp=100;

        // Maximum Pressure for which database is available
        static const int maxPress=9890000.00;

       //Number of temperature steps
        static const int nts=(maxTemp-minTemp)/stepTemp;




public:

    //- Runtime type information
    TypeName("ngrAbsorptionEmission");


    // Constructors

        //- Construct from components
        ngrAbsorptionEmission
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~ngrAbsorptionEmission();


    // Member Functions

        // Access

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
			
        void correct
        (
            volScalarField& a,
            PtrList<volScalarField>& aLambda
        ) const;

    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 // End namespace absorptionEmissionModels
} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
