/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "oqradiativeIntensityRay.H"
#include "fvm.H"
#include "oqfvDOM.H"
#include "constants.H"

using namespace Foam::constant;

const Foam::word
Foam::radiation::oqradiativeIntensityRay::intensityPrefix("ILambda");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::oqradiativeIntensityRay::oqradiativeIntensityRay
(
    const oqfvDOM& dom,
    const fvMesh& mesh,
    const scalar phi,
    const scalar theta,
    const scalar deltaPhi,
    const scalar deltaTheta,
    const label nLambda,
    const absorptionEmissionModel& absorptionEmission,
    const blackBodyEmission& blackBody,
    const label rayId
)
:
    dom_(dom),
    mesh_(mesh),
    absorptionEmission_(absorptionEmission),
    blackBody_(blackBody),
    I_
    (
        IOobject
        (
            "I" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("I", dimMass/pow3(dimTime), 0.0)
    ),
    qr_
    (
        IOobject
        (
            "qr" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("qr", dimMass/pow3(dimTime), 0.0)
    ),
    qin_
    (
        IOobject
        (
            "qin" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("qin", dimMass/pow3(dimTime), 0.0)
    ),
    qem_
    (
        IOobject
        (
            "qem" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("qem", dimMass/pow3(dimTime), 0.0)
    ),
    d_(Zero),
    dAve_(Zero),
    theta_(theta),
    phi_(phi),
    omega_(0.0),
    nLambda_(nLambda),
    ILambda_(nLambda),
    myRayId_(rayId)
{
    scalar sinTheta = Foam::sin(theta);
    scalar cosTheta = Foam::cos(theta);
    scalar sinPhi = Foam::sin(phi);
    scalar cosPhi = Foam::cos(phi);

    omega_ = 2.0*sinTheta*Foam::sin(deltaTheta/2.0)*deltaPhi;
    d_ = vector(sinTheta*sinPhi, sinTheta*cosPhi, cosTheta);
    
    Info<<"the d_ are "<<d_<<endl;
    
    dAve_ = vector
    (
        sinPhi
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*theta)
       *Foam::sin(deltaTheta)),
        cosPhi
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*theta)
       *Foam::sin(deltaTheta)),
        0.5*deltaPhi*Foam::sin(2.0*theta)*Foam::sin(deltaTheta)
    );
     
   //  Info<<"the dAve_ are "<<dAve_<<endl;

    if (mesh_.nSolutionD() == 2)
    {
        vector meshDir(Zero);
        if (dom_.meshOrientation() != vector::zero)
        {
            meshDir = dom_.meshOrientation();
        }
        else
        {
            for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
            {
                if (mesh_.geometricD()[cmpt] == -1)
                {
                    meshDir[cmpt] = 1;
                }
            }
        }
        const vector normal(vector(0, 0, 1));

        const tensor coordRot = rotationTensor(normal, meshDir);

        dAve_ = coordRot & dAve_;
        d_ = coordRot & d_;
    }
    else if (mesh_.nSolutionD() == 1)
    {
        vector meshDir(Zero);
        if (dom_.meshOrientation() != vector::zero)
        {
            meshDir = dom_.meshOrientation();
        }
        else
        {
            for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
            {
                if (mesh_.geometricD()[cmpt] == 1)
                {
                    meshDir[cmpt] = 1;
                }
            }
        }
        const vector normal(vector(1, 0, 0));

        dAve_ = (dAve_ & normal)*meshDir;
        d_ = (d_ & normal)*meshDir;
    }

    autoPtr<volScalarField> IDefaultPtr;

    forAll(ILambda_, lambdaI)
    {
        IOobject IHeader
        (
            intensityPrefix + "_" + name(rayId) + "_" + name(lambdaI),
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        // Check if field exists and can be read
        if (IHeader.typeHeaderOk<volScalarField>(true))
        {
            ILambda_.set
            (
                lambdaI,
                new volScalarField(IHeader, mesh_)
            );
        }
        else
        {
            // Demand driven load the IDefault field
            if (!IDefaultPtr.valid())
            {
                IDefaultPtr.reset
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "IDefault",
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh_
                    )
                );
            }

            // Reset the MUST_READ flag
            IOobject noReadHeader(IHeader);
            noReadHeader.readOpt() = IOobject::NO_READ;

            ILambda_.set
            (
                lambdaI,
                new volScalarField(noReadHeader, IDefaultPtr())
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::oqradiativeIntensityRay::~oqradiativeIntensityRay()
{}

// Info<< "End\n" << endl;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::radiation::oqradiativeIntensityRay::correct()
{
    // Reset boundary heat flux to zero
    qr_.boundaryFieldRef() = 0.0;

    scalar maxResidual = -GREAT;
  
   // scalarField f1Cells =I_.internalField();
      
    forAll(ILambda_, lambdaI)
    {
       // forAll(f1Cells,celli)
       
        const volScalarField& k = dom_.aLambda(lambdaI);

        const surfaceScalarField Ji(dAve_ & mesh_.Sf());
               
        fvScalarMatrix IiEq
        (
            fvm::div(Ji, ILambda_[lambdaI], "div(Ji,Ii_h)")
            
              + 
              
           fvm::Sp(k*omega_, ILambda_[lambdaI])
        ==
            1.0/constant::mathematical::pi*omega_
           *(
                (k - absorptionEmission_.aDisp(lambdaI))
               *blackBody_.bLambda(lambdaI)

              + absorptionEmission_.E(lambdaI)/4
            )
        );

        
      IiEq.relax();
       // Info<<"the intensity values are"<< ILambda_[lambdaI]<<endl;

        const solverPerformance ILambdaSol = solve
        (
            IiEq,
            mesh_.solver("Ii")
        );
        
    const fvPatchList& patches = mesh_.boundary();
     
    forAll(ILambda_, lambdaI)
    { 
       forAll(ILambda_[lambdaI].internalField(),celli)
        {
             scalar I=ILambda_[lambdaI].internalField()[celli];
            // Info<<"i am in for loop"<<endl;   
          //if(f1Cells[celli] < 0.0)
          if(I < 1e0)
          {  
             // Info<<"I am in if loop"<<endl;
              ILambda_[lambdaI][celli]=0.0;
           }
        }
        
        forAll(patches, patchi)
		{
			//const fvPatch& currPatch = patches[patchi];
			forAll(ILambda_[lambdaI].boundaryField(), patchi)
			{
       scalarField Int_suf = ILambda_[lambdaI].boundaryField()[patchi];//collecting the intensity for the patch
				forAll(Int_suf, facei)
				{
					scalar Int = Int_suf[facei]; 
					/*if(myRayId_ == 6)
					{
						Info<<"BoundaryField::"<<endl;
						Info<<"s_I = "<<s_I<<endl;
					}*/
					
					if(Int < 0.0)
					{
						ILambda_[lambdaI].boundaryFieldRef()[patchi][facei] = 0.0;
					}
				}
			}			
		}
  	
    }  //end of the for loop 
    
    
        const scalar initialRes = ILambdaSol.initialResidual()*omega_/dom_.omegaMax();

        maxResidual = max(initialRes, maxResidual);
    }
  
 
    return maxResidual;
}


void Foam::radiation::oqradiativeIntensityRay::addIntensity()
{
    I_ = dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);

    forAll(ILambda_, lambdaI)
    {
        
        I_ += ILambda_[lambdaI];
    }
}


// ************************************************************************* //
