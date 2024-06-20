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
    Foam::radiation::nonGrayAbsorptionEmission

Group
    grpRadiationAbsorptionEmissionSubModels

Description

    nonGrayAbsorptionEmission radiation absorption/emission for solids.


SourceFiles
    nonGrayAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef nonGrayAbsorptionEmission_H
#define nonGrayAbsorptionEmission_H
#include "absorptionCoeffs.H"
#include "absorptionEmissionModel.H"
#include "fluidThermo.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class nonGrayAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class nonGrayAbsorptionEmission
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

       //- Thermo package
      //  const fluidThermo& thermo_;

      // Total set of moleFractions available in dataBase
      static const int nMoleFractions=12;
			
      // Minimum Temperature at which database is available
      static const int minTemp=300;

      // Temperature step for which database is available
       static const int stepTemp=100;



public:

    //- Runtime type information
    TypeName("nonGrayAbsorptionEmission");


    // Constructors

        //- Construct from components
        nonGrayAbsorptionEmission
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~nonGrayAbsorptionEmission();


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
