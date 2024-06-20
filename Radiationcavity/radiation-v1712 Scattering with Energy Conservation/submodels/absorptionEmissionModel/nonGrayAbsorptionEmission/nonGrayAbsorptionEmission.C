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

#include "nonGrayAbsorptionEmission.H"
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
    	
	defineTypeNameAndDebug(nonGrayAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            nonGrayAbsorptionEmission,
            dictionary
        );
    
     }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::nonGrayAbsorptionEmission::
nonGrayAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    absCoeffs_(maxBands_),
    emiCoeffs_(maxBands_),
    //thermo_(mesh.lookupObject<fluidThermo>(basicThermo::dictName)), 
    nBands_(32)
{
    coeffsDict_.lookup("absorptivity") >> absCoeffs_;
    coeffsDict_.lookup("emissivity") >> emiCoeffs_;
    nBands_ = absCoeffs_.size();
    coeffsDict_.lookup("weight")>>weights_;
    nWeights_ =weights_.size();
		
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::nonGrayAbsorptionEmission::
~nonGrayAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::nonGrayAbsorptionEmission::aCont
(
    const label maxBands_
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
	
//const basicSpecieMixture& mixture = dynamic_cast<const basicSpecieMixture&>(thermo_);

  // const volScalarField& T = thermo_.T();
   // const volScalarField& p = thermo_.p();
		
         
	const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
	const volScalarField& H2O = mesh_.lookupObject<volScalarField>("H2O");
	int iter=0;int iterprev=0; int iterpost=1;
	forAll(mesh_.V(),celli)
   	{

	  volScalarField itprev
   	(
        IOobject
        (
            "itprev",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("itprev",dimless, ROOTVSMALL)
    );

		volScalarField itpost
   	(
        IOobject
        (
            "itpost",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("itpost",dimless, ROOTVSMALL)
    );

     //  scalar moleFractions[12]={0.01,0.02,0.03,0.04,0.05,0.1,0.15,0.2,0.25,0.5,0.75,1};
			 scalarField intermole(12);	
			 intermole[0]=0.01;
	     intermole[1]=0.02;
			 intermole[2]=0.03;
			 intermole[3]=0.04;
			 intermole[4]=0.05;
			 intermole[5]=0.1;
			 intermole[6]=0.15;
			 intermole[7]=0.2;
			 intermole[8]=0.25;
		   intermole[9]=0.5;
			 intermole[10]=0.75;
			 intermole[11]=1;
/*	int m=0;
	while((H2O[celli]-moleFractions[m])>=0)
	{
	     m++;
	}
	m--;*/
    linearInterpolationWeights interpolator
    //splineInterpolationWeights interpolator
    (
        intermole
    );
    labelList indices;

    scalarField weights;

    interpolator.valueWeights(H2O[celli], indices, weights);

	 iter=1; iterprev=1; iterpost=1;

   iter=(T[celli]-minTemp)/stepTemp;

	 itprev[celli]=(iter*nMoleFractions+indices[0])*nWeights_+maxBands_;

	 itpost[celli]=(iter*nMoleFractions+indices[1])*nWeights_+maxBands_;

	 iterprev=itprev[celli];

	 iterpost=itpost[celli];

	 a[celli]=absCoeffs_[iterprev]+(1-weights[0])*(absCoeffs_[iterpost]-absCoeffs_[iterprev]);

	/* iter=(T[celli]-300)/100;
	 it[celli]=(iter*12+m)*nWeights_+maxBands_;
	 iter1=it[celli];
	 a[celli]=absCoeffs_[iter1];*/

	}

    return ta;
}



Foam::tmp<Foam::volScalarField>
Foam::radiation::nonGrayAbsorptionEmission::eCont
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
Foam::radiation::nonGrayAbsorptionEmission::ECont
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

void Foam::radiation::nonGrayAbsorptionEmission::correct
(
    volScalarField& a,
    PtrList<volScalarField>& aLambda
) const

{
   a = dimensionedScalar("zero", dimless/dimLength, 0.0);

		
  for (label j=0; j<nWeights_; j++)
    {
	
	
        aLambda[j].primitiveFieldRef() = this->a(j);

	volScalarField bLambda
   	(
        IOobject
        (
            "bLambda",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
       this->a(j)
       );
		
aLambda[j].boundaryFieldRef()=bLambda.boundaryFieldRef();

a.primitiveFieldRef() = aLambda[j].primitiveField();
			
a.boundaryFieldRef()=bLambda.boundaryFieldRef();
		
	}	

}


// ************************************************************************* //
