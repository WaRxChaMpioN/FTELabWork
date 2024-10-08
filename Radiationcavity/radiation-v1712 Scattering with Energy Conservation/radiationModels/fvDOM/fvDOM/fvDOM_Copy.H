/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::radiation::fvDOM

Group
    grpRadiationModels

Description

    Finite Volume Discrete Ordinates Method. Solves the RTE equation for n
    directions in a participating media, not including scatter.

    Available absorption models:
    - constantAbsorptionEmission
    - greyMeanAbsoprtionEmission
    - wideBandAbsorptionEmission

Usage
    \verbatim
        fvDOMCoeffs
        {
            nPhi        4;          // azimuthal angles in PI/2 on X-Y.
                                    //(from Y to X)
            nTheta      0;          // polar angles in PI (from Z to X-Y plane)
            convergence 1e-3;       // convergence criteria for radiation
                                    //iteration
            maxIter     4;          // maximum number of iterations
             meshOrientation    (1 1 1); //Mesh ortientation used for 2D and 1D
        }

        solverFreq   1; // Number of flow iterations per radiation iteration
    \endverbatim

    The total number of solid angles is  4*nPhi*nTheta in 3-D.

    Operating modes:
    - 1-D:
      - ray directions are on X, Y or Z
      - \c nPhi and \c nTheta entries are ignored
      - \c meshOrientation vector can be used for any other 1-D direction.
    - 2-D:
      - ray directions are on X-Y, X-Z or Y-Z planes
      - only the \c nPhi entry is considered
      - \c meshOrientation vector can be used for non-aligned planes
        specifying the plane normal vector.
    - 3-D:
      - rays geberated in 3-D using the \c nPhi and \c nTheta entries
      - \c meshOrientation vector is not applicable.

SourceFiles
    fvDOM.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_fvDOM_H
#define radiation_fvDOM_H

#include "radiativeIntensityRay.H"
#include "radiationModel.H"
#include "fvMatrices.H"
#include "solarLoad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
// Forward declaration of classes
class scatterModel;
class radiativeIntensityRay;

/*---------------------------------------------------------------------------*\
                           Class fvDOM Declaration
\*---------------------------------------------------------------------------*/

class fvDOM
:
    public radiationModel
{
    // Private data


        //- Incident radiation  [W/m2]
        volScalarField G_;

  	//- Incident radiation  [W/m2]
        volScalarField Gi;
       
       //Sum the intensities scatter form the direction other than i
       volScalarField I_Sum;
        
	volScalarField delq;

        //- Total radiative heat flux [W/m2]
        volScalarField qr_;

         //- Emmited radiative heat flux [W/m2]
        volScalarField qem_;

	//- Emmited radiative heat flux [W/m2]
        volScalarField qRu;

	//- Emmited radiative heat flux [W/m2]
        volScalarField qRP;

        //- Incidet radiative heat flux [W/m2]
        volScalarField qin_;

        //- Total absorption coefficient [1/m]
        volScalarField a_;

        //- Number of solid angles in theta
        label nTheta_;  //label is int 32 or int 64

        //- Number of solid angles in phi
        label nPhi_ ;

        //- Total number of rays (1 per direction)
        label nRay_;

        //- Number of wavelength bands
        label nLambda_;

        //- Wavelength total absorption coefficient [1/m]
        PtrList<volScalarField> aLambda_; //array object alambda
        
        //- Black body
        blackBodyEmission blackBody_;

        //- List of pointers to radiative intensity rays
        PtrList<radiativeIntensityRay> IRay_;  //array object of class radiativeintensity

        //- Convergence criterion
        scalar convergence_;

        //- Maximum number of iterations
        label maxIter_;

        //- Maximum omega weight
        scalar omegaMax_;

        //- Use Solar Load model
        bool useSolarLoad_;

        //- Solar load radiation model
        autoPtr<solarLoad> solarLoad_;

        //- Mesh orientation vector
        vector meshOrientation_;


    // Private Member Functions

        //- Initialise
        void initialise();

        //- No copy Constructor
        fvDOM(const fvDOM&);

        //- No copy assignment
        void operator=(const fvDOM&);

        //- Update black body emission
        void updateBlackBodyEmission();


public:

    //- Runtime type information
    TypeName("fvDOM");


    // Constructors

        //- Construct from components
        fvDOM(const volScalarField& T);

        //- Construct from components
        fvDOM(const dictionary& dict, const volScalarField& T);


    //- Destructor
    virtual ~fvDOM();


    // Member functions

        // Edit

            //- Solve radiation equation(s)
            void calculate();
            
            //void sum_Zero();

            //- Read radiation properties dictionary
            bool read();

            //- Update G and calculate total heat flux on boundary
            void updateG();

            //- Set the rayId and lambdaId from by decomposing an intensity
            //  field name
            void setRayIdLambdaId
            (
                const word& name,
                label& rayId,
                label& lambdaId
            ) const;

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const;
					
	    void delqr();



        // Access

            //- Ray intensity for rayI
            inline const radiativeIntensityRay& IRay(const label rayI) const;

            //- Ray intensity for rayI and lambda bandwidth
            inline const volScalarField& IRayLambda
            (
                const label rayI,
                const label lambdaI
            ) const;

            //- Number of angles in theta
            inline label nTheta() const;

            //- Number of angles in phi
            inline label nPhi() const;

            //- Number of rays
            inline label nRay() const;

            //- Number of wavelengths
            inline label nLambda() const;
            
            //return the Sum of the intensities scatter from the other direction into the ray i direction 
            inline const volScalarField& return_ScatterSum() const;

            //- Const access to total absorption coefficient
            inline const volScalarField& a() const;

            //- Const access to wavelength total absorption coefficient
            inline const volScalarField& aLambda(const label lambdaI) const;

            //- Const access to incident radiation field
            inline const volScalarField& G() const;

            //- Const access to total radiative heat flux field
            inline const volScalarField& qr() const;

            //- Const access to incident radiative heat flux field
            inline const volScalarField& qin() const;

            //- Const access to emitted radiative heat flux field
            inline const volScalarField& qem() const;

            //- Const access to black body
            inline const blackBodyEmission& blackBody() const;

            //- Return omegaMax
            inline scalar omegaMax() const;

            //- Return meshOrientation
            inline vector meshOrientation() const;

	    //- Return divergence of radiative heat flux
	    //	inline const volScalarField& delqr() const;
	    //- Number of bands
            label nBands() const;
					
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvDOMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
