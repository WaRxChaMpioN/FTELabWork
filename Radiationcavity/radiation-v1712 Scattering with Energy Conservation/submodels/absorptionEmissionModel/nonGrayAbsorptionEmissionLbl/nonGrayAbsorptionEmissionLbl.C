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

\*---------------------------------------------------------------------------*/

#include "nonGrayAbsorptionEmissionLbl.H"
#include "addToRunTimeSelectionTable.H"
#include "basicSpecieMixture.H"
#include "HashTable.H"
#include "scalarIndList.H"
#include "IOdictionary.H"
#include "linearInterpolationWeights.H"
#include "splineInterpolationWeights.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
    	
	defineTypeNameAndDebug(nonGrayAbsorptionEmissionLbl, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            nonGrayAbsorptionEmissionLbl,
            dictionary
        );
    
     }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::nonGrayAbsorptionEmissionLbl::
nonGrayAbsorptionEmissionLbl
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    absCoeffs_(maxBands_),
    emiCoeffs_(maxBands_),
    nBands_(maxBands_)
{
    coeffsDict_.lookup("absorptivity") >> absCoeffs_;  //absCoeffs is a scalar list so it store aal the value as a array
    coeffsDict_.lookup("emissivity") >> emiCoeffs_; //emiCoeffs is a scalar list so it store aal the value as a array
    nBands_ = absCoeffs_.size(); //Deifne the size for the nBands according the value in absCoeff file
    /*coeffsDict_.lookup("weight")>>weights_; 
    nWeights_ =weights_.size(); //Deifne the size for the nBands according the value in absCoeff file*/
		
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::nonGrayAbsorptionEmissionLbl::
~nonGrayAbsorptionEmissionLbl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::nonGrayAbsorptionEmissionLbl::aCont
(
    const label bandI
) const
{
	
  	tmp<volScalarField> ta
        (
        new volScalarField
        (
            IOobject
            (
                "a",
               mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a",dimless/dimLength,0.0)
					
        )
        );
	
	scalarField& a = ta.ref().primitiveFieldRef();
	
	forAll(mesh_.V(),celli)
   	{
   	    a = absCoeffs_[bandI];

	}

    return ta;
}



Foam::tmp<Foam::volScalarField>
Foam::radiation::nonGrayAbsorptionEmissionLbl::eCont
(
    const label bandI
) const
{
    tmp<volScalarField> te
    (
        new volScalarField
        (
            IOobject
            (
                "e",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("e", dimless/dimLength, emiCoeffs_[bandI])
        )
    );

    return te;

}


Foam::tmp<Foam::volScalarField>
Foam::radiation::nonGrayAbsorptionEmissionLbl::ECont
(
    const label bandI
) const
{
    tmp<volScalarField> E
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
           dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    return E;
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::nonGrayAbsorptionEmissionLbl::a(const label bandI) const
{
    return  aCont(bandI);
}

void Foam::radiation::nonGrayAbsorptionEmissionLbl::correct
(
    volScalarField& a,
    PtrList<volScalarField>& aj
) const

{		
    for (label j=0; j<nBands_; j++)
    {
        a = this->a(j);
        aj[j] = a;
    }
		
}



// ************************************************************************* //
