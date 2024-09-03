/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "oqfvDOM.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "constants.H"
#include "unitConversion.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(oqfvDOM, 0);
        addToRadiationRunTimeSelectionTables(oqfvDOM);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::oqfvDOM::rotateInitialRays(const vector& sunDir)
{
    // Rotate Y spherical cordinates to Sun direction.
    // Solid angles on the equator are better fit for planar radiation
    const tensor coordRot = rotationTensor(vector(0, 1, 0), sunDir);

    forAll(IRay_, rayId)
    {
        IRay_[rayId].dAve() = coordRot & IRay_[rayId].dAve();
        IRay_[rayId].d() = coordRot & IRay_[rayId].d();
    }
}


void Foam::radiation::oqfvDOM:: alignClosestRayToSun(const vector& sunDir)
{
    label SunRayId(-1);
    scalar maxSunRay = -GREAT;

    // Looking for the ray closest to the Sun direction
    forAll(IRay_, rayId)
    {
        const vector& iD = IRay_[rayId].d();
        scalar dir = sunDir & iD;
        if (dir > maxSunRay)
        {
            maxSunRay = dir;
            SunRayId = rayId;
        }
    }

    // Second rotation to align colimated radiation with the closest ray
    const tensor coordRot = rotationTensor(IRay_[SunRayId].d(), sunDir);

    forAll(IRay_, rayId)
    {
        IRay_[rayId].dAve() = coordRot & IRay_[rayId].dAve();
        IRay_[rayId].d() = coordRot & IRay_[rayId].d();
    }

    Info << "Sun direction : " << sunDir << nl << endl;
    Info << "Sun ray ID : " << SunRayId << nl << endl;
}


void Foam::radiation::oqfvDOM::updateRaysDir()
{
    solarCalculator_->correctSunDirection();
    const vector sunDir = solarCalculator_->direction();

    // First iteration
    if (updateTimeIndex_ == 0)
    {
        rotateInitialRays(sunDir);
        alignClosestRayToSun(sunDir);
    }
    else if (updateTimeIndex_ > 0)
    {
        alignClosestRayToSun(sunDir);
    }
}


void Foam::radiation::oqfvDOM::initialise()
{
    coeffs_.readIfPresent("useExternalBeam", useExternalBeam_);

    if (useExternalBeam_)
    {
        spectralDistributions_.reset
        (
            Function1<scalarField>::New
            (
                "spectralDistribution",
                coeffs_,
                &mesh_
            )
        );

        spectralDistribution_ =
            spectralDistributions_->value(mesh_.time().timeOutputValue());

        spectralDistribution_ =
            spectralDistribution_/sum(spectralDistribution_);

        const dictionary& solarDict = this->subDict("solarCalculatorCoeffs");
        solarCalculator_.reset(new solarCalculator(solarDict, mesh_));

        if (mesh_.nSolutionD() != 3)
        {
            FatalErrorInFunction
                << "External beam model only available in 3D meshes "
                << abort(FatalError);
        }

        if (solarCalculator_->diffuseSolarRad() > 0)
        {
            FatalErrorInFunction
                << "External beam model does not support Diffuse "
                << "Solar Radiation. Set diffuseSolarRad to zero"
                << abort(FatalError);
        }
        if (spectralDistribution_.size() != nLambda_)
        {
            FatalErrorInFunction
                << "The epectral energy distribution has different bands "
                << "than the absoprtivity model "
                << abort(FatalError);
        }
    }

    // 3D
    if (mesh_.nSolutionD() == 3)
    {
        nRay_ = 4*nPhi_*nTheta_;

        IRay_.setSize(nRay_);

        const scalar deltaPhi = pi/(2*nPhi_);
        const scalar deltaTheta = pi/nTheta_;

        label i = 0;

        for (label n = 1; n <= nTheta_; n++)
        {
            for (label m = 1; m <= 4*nPhi_; m++)
            {
                scalar thetai = (2*n - 1)*deltaTheta/2.0;
                scalar phii = (2*m - 1)*deltaPhi/2.0;

                IRay_.set
                (
                    i,
                    new oqradiativeIntensityRay
                    (
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        *absorptionEmission_,
                        blackBody_,
                        i
                    )
                );
                i++;
            }
        }
    }
    // 2D
    else if (mesh_.nSolutionD() == 2)
    {
        const scalar thetai = piByTwo;
        const scalar deltaTheta = pi;
        nRay_ = 4*nPhi_;
        IRay_.setSize(nRay_);
        const scalar deltaPhi = pi/(2.0*nPhi_);
        label i = 0;
        for (label m = 1; m <= 4*nPhi_; m++)
        {
            const scalar phii = (2*m - 1)*deltaPhi/2.0;
            IRay_.set
            (
                i,
                new oqradiativeIntensityRay
                (
                    *this,
                    mesh_,
                    phii,
                    thetai,
                    deltaPhi,
                    deltaTheta,
                    nLambda_,
                    *absorptionEmission_,
                    blackBody_,
                    i
                )
            );
            i++;
        }
    }
    // 1D
    else
    {
        const scalar thetai = piByTwo;
        const scalar deltaTheta = pi;
        nRay_ = 2;
        IRay_.setSize(nRay_);
        const scalar deltaPhi = pi;
        label i = 0;
        for (label m = 1; m <= 2; m++)
        {
            const scalar phii = (2*m - 1)*deltaPhi/2.0;
            IRay_.set
            (
                i,
                new oqradiativeIntensityRay
                (
                    *this,
                    mesh_,
                    phii,
                    thetai,
                    deltaPhi,
                    deltaTheta,
                    nLambda_,
                    *absorptionEmission_,
                    blackBody_,
                    i
                )
            );
            i++;
        }
    }


    // Construct absorption field for each wavelength
    forAll(aLambda_, lambdaI)
    {
        aLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "aLambda_" + Foam::name(lambdaI) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                a_
            )
        );
    }

    Info<< "oqfvDOM : Allocated " << IRay_.size()
        << " rays with average orientation:" << nl;

    if (useExternalBeam_)
    {
        // Rotate rays for Sun direction
        updateRaysDir();
    }

    scalar totalOmega = 0;
    forAll(IRay_, rayId)
    {
        if (omegaMax_ <  IRay_[rayId].omega())
        {
            omegaMax_ = IRay_[rayId].omega();
        }
        totalOmega += IRay_[rayId].omega();
        Info<< '\t' << IRay_[rayId].I().name() << " : " << "dAve : "
            << '\t' << IRay_[rayId].dAve() << " : " << "omega : "
            << '\t' << IRay_[rayId].omega() << " : " << "d : "
            << '\t' << IRay_[rayId].d() << nl;
    }

    Info << "Total omega : " << totalOmega << endl;

    Info<< endl;

    coeffs_.readIfPresent("useSolarLoad", useSolarLoad_);

    if (useSolarLoad_)
    {
        // if (useExternalBeam_)
        // {
        //     FatalErrorInFunction
        //         << "External beam with oqfvDOM can not be used "
        //         << "with the solar load model"
        //         << abort(FatalError);
        // }
        // const dictionary& solarDict = this->subDict("solarLoadCoeffs");
        // solarLoad_.reset(new solarLoad(solarDict, T_));

        // if (solarLoad_->nBands() != this->nBands())
        // {
        //     FatalErrorInFunction
        //         << "Requested solar radiation with oqfvDOM. Using "
        //         << "different number of bands for the solar load is not allowed"
        //         << abort(FatalError);
        // }

        // Info<< "Creating Solar Load Model " << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::oqfvDOM::oqfvDOM(const volScalarField& T)
:
    radiationModel(typeName, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    Gi                          //aded
    (
        IOobject
        (
            "Gi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Gi", dimensionSet(1, -1, -3, 0, 0), 0.0)
    ),
    delqr_                          //aded
    (
        IOobject
        (
            "delqr_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("delq", dimensionSet(1, -1, -3, 0, 0), 0.0)
    ),
		qRu                 //aded
    (
        IOobject
        (
            "qRu",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qRu", dimensionSet(1, -1, -3, 0, 0), 0.0)
    ),
    
    I_Sum                   //aded
    (
        IOobject
        (
            "I_Sum",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("I_Sum", dimMass/pow3(dimTime), 0.0)
    ),
	
    qRP                     //aded
    (
        IOobject
        (
            "qRP",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
       dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0)
    ),
    
    qrB                     //aded
    (
        IOobject
        (
            "qrB",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
       dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0)
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qem_
    (
        IOobject
        (
            "qem",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qin_
    (
        IOobject
        (
            "qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),

    nTheta_(coeffs_.get<label>("nTheta")),
    nPhi_(coeffs_.get<label>("nPhi")),
    nRay_(0),
    nLambda_(readLabel(coeffs_.lookup("nWeights"))), //Uncomment it when not using the LBL case   //aded
    // nLambda_(absorptionEmission_->nBands()), // Uncomment it when using the LBL case   //removed
    aLambda_(nLambda_),
    blackBody_(nLambda_, T),
    IRay_(0),
    tolerance_
    (
        coeffs_.getOrDefaultCompat<scalar>
        (
            "tolerance",
            {{"convergence", 1712}},
            0
        )
    ),
    maxIter_(coeffs_.getOrDefault<label>("maxIter", 50)),
    omegaMax_(0),
    useSolarLoad_(false),
    solarLoad_(),
    meshOrientation_
    (
        coeffs_.getOrDefault<vector>("meshOrientation", Zero)
    ),
    useExternalBeam_(false),
    spectralDistribution_(),
    spectralDistributions_(),
    solarCalculator_(),
    updateTimeIndex_(0)
{
    initialise();
}


Foam::radiation::oqfvDOM::oqfvDOM
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    Gi              //aded
    (
        IOobject
        (
            "Gi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Gi", dimensionSet(1, -1, -3, 0, 0), 0.0)
    ),
    delqr_                  //aded
    (
        IOobject
        (
            "delqr_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("delqr_", dimensionSet(1, -1, -3, 0, 0), 0.0)
    ),
    qRu                 //aded
    (
        IOobject
        (
            "qRu",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qRu", dimensionSet(1, -1, -3, 0, 0), 0.0)
    ),
    
    qrB                 //aded
    (
        IOobject
        (
            "qrB",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
       dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0)
    ),
    I_Sum                   //aded
    (
        IOobject
        (
            "I_Sum",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("I_Sum", dimMass/pow3(dimTime), 0.0)
    ),
    
    qRP                 //aded
    (
        IOobject
        (
            "qRP",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0)
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qem_
    (
        IOobject
        (
            "qem",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qin_
    (
        IOobject
        (
            "qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    nTheta_(coeffs_.get<label>("nTheta")),
    nPhi_(coeffs_.get<label>("nPhi")),
    nRay_(0),
    nLambda_(readLabel(coeffs_.lookup("nWeights"))),   //Uncomment it when not using the LBL case  //aded
    // nLambda_(absorptionEmission_->nBands()),    // Uncomment it when using the LBL case   //commented
    aLambda_(nLambda_),
    blackBody_(nLambda_, T),
    IRay_(0),
    tolerance_
    (
        coeffs_.getOrDefaultCompat<scalar>
        (
            "tolerance",
            {{"convergence", 1712}},
            0
        )
    ),
    maxIter_(coeffs_.getOrDefault<label>("maxIter", 50)),
    omegaMax_(0),
    useSolarLoad_(false),
    solarLoad_(),
    meshOrientation_
    (
        coeffs_.getOrDefault<vector>("meshOrientation", Zero)
    ),
    useExternalBeam_(false),
    spectralDistribution_(),
    spectralDistributions_(),
    solarCalculator_(),
    updateTimeIndex_(0)
{
    initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::oqfvDOM::read()
{
    if (radiationModel::read())
    {
        // Only reading solution parameters - not changing ray geometry
        coeffs_.readIfPresentCompat
        (
            "tolerance", {{"convergence", 1712}}, tolerance_
        );
        coeffs_.readIfPresent("maxIter", maxIter_);

        return true;
    }

    return false;
}


void Foam::radiation::oqfvDOM::calculate()
{
    absorptionEmission_->correct(a_, aLambda_);

    updateBlackBodyEmission();

    if (useSolarLoad_)
    {
        solarLoad_->calculate();
    }

    if (useExternalBeam_)
    {
        // switch (solarCalculator_->sunDirectionModel())
        // {
        //     case solarCalculator::mSunDirConstant:
        //     {
        //         break;
        //     }
        //     case solarCalculator::mSunDirTracking:
        //     {
        //         label updateIndex = label
        //         (
        //             mesh_.time().value()
        //            /solarCalculator_->sunTrackingUpdateInterval()
        //         );

        //         if (updateIndex > updateTimeIndex_)
        //         {
        //             Info << "Updating Sun position..." << endl;
        //             updateTimeIndex_ = updateIndex;
        //             updateRaysDir();
        //         }
        //         break;
        //     }
        // }
    }

    // Set rays convergence false
    List<bool> rayIdConv(nRay_, false);

    scalar maxResidual = 0;
    label radIter = 0;
    do
    {
        Info<< "Radiation solver iter: " << radIter << endl;

        radIter++;
        maxResidual = 0;

        //Update I
        forAll(IRay_,rayI)              //aded
        {
            IRay_[rayI].updateI_Old();
        }

        forAll(IRay_, rayI)
        {

            I_Sum = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);   //aded
            
            label new_rayI = rayI;  //aded

            forAll(IRay_,rayI)          //aded
            {
                if(rayI != new_rayI)
                {
                    I_Sum += IRay_[rayI].returnI_Old();
                }
            }

            if (!rayIdConv[rayI])
            {
                scalar maxBandResidual = IRay_[rayI].correct();
                maxResidual = max(maxBandResidual, maxResidual);

                if (maxBandResidual < tolerance_)
                {
                    rayIdConv[rayI] = true;
                }
            }
        }

    } while (maxResidual > tolerance_ && radIter < maxIter_);

    updateG();
    qr_Domain();           //aded 
    qr_Boundary();      //aded
    energyCheck();      //aded
}


Foam::tmp<Foam::volScalarField> Foam::radiation::oqfvDOM::Rp() const
{

  List<scalar> aj (coeffs_.lookup("weight"));   //aded

    // Construct using contribution from first frequency band
    auto tRp = volScalarField::New
    (
        "Rp",
        IOobject::NO_REGISTER,
        (
            4
          * physicoChemical::sigma
          * (aLambda_[0] - absorptionEmission_->aDisp(0)())
          * blackBody_.deltaLambdaT(T_, absorptionEmission_->bands(0))
        )
    );
    auto& Rp = tRp.ref();

    // Add contributions over remaining frequency bands
    for (label j=1; j < nLambda_; j++)
    {
        Rp +=
        (
            4
           *physicoChemical::sigma
           *(aLambda_[j] - absorptionEmission_->aDisp(j)())
           *blackBody_.deltaLambdaT(T_, absorptionEmission_->bands(j))
        );
    }

    return tRp;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiation::oqfvDOM::Ru() const
{
    auto tRu = DimensionedField<scalar, volMesh>::New
    (
        "Ru",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimensionSet(1, -1, -3, 0, 0), Zero)
    );
    auto& Ru = tRu.ref();

    // Sum contributions over all frequency bands
    for (label j=0; j < nLambda_; j++)
    {
        // Compute total incident radiation within frequency band
        tmp<DimensionedField<scalar, volMesh>> Gj
        (
            IRay_[0].ILambda(j)()*IRay_[0].omega()
        );

        for (label rayI=1; rayI < nRay_; rayI++)
        {
            Gj.ref() += IRay_[rayI].ILambda(j)()*IRay_[rayI].omega();
        }

        Ru += (aLambda_[j]() - absorptionEmission_->aDisp(j)()())*Gj
             - absorptionEmission_->ECont(j)()();
    }

    return tRu;
}


void Foam::radiation::oqfvDOM::updateBlackBodyEmission()
{
    for (label j=0; j < nLambda_; j++)
    {
        blackBody_.correct(j, absorptionEmission_->bands(j));
    }
}

void Foam::radiation::oqfvDOM::qr_Domain()            //aded
{
   qr_Domain_ = 0.0; //- need to initialise becasue it take previous iteration value
   
   forAll(mesh_.cells(),celli)  //for(label celli = 0; celli<mesh_.nCells(); celli++)  
   {
       qr_Domain_ += delqr_[celli]*mesh_.V()[celli];
   } 
}

void Foam::radiation::oqfvDOM::qr_Boundary()          //aded
{  
   qr_Boundary_ = 0.0; //- need to initialise becasue it take previous iteration value
   forAll(mesh_.boundary(), patchi)
   {
      scalar qrface_ = 0.0;
      const scalarField& qrpatch = qr_.boundaryField()[patchi];
      const scalarField& area = mesh_.boundary()[patchi].magSf();   
      if(qr_.boundaryField()[patchi].size()>0)
      {
          forAll(mesh_.boundary()[patchi],facei)
          {
                qrface_ += qrpatch[facei]*area[facei];
          }
      } 
      qr_Boundary_ += qrface_; 
   } 
}
void Foam::radiation::oqfvDOM::energyCheck()              //aded
{
    Info<<"qr_Boundary = "<<qr_Boundary_<<endl;
    Info<<"qr_Domain = "<<qr_Domain_<<endl;
    
    scalar Diff = abs(qr_Domain_ - qr_Boundary_);

    if(abs(qr_Domain_)>1)
    {
        scalar error = abs((Diff/qr_Domain_)*100);
        Info<<"Error = "<<error<<endl;
        if(error < 0.1)
        {
            Info<<"Energy is conserved"<<endl;
        }
    }
    
    else if(Diff < 20.0)
    {
        Info<<"Energy is conserved"<<endl;
    }
}
void Foam::radiation::oqfvDOM::qr_Cell()                  //aded
{
   forAll(mesh_.cells(),celli)  //for(label celli = 0; celli<mesh_.nCells(); celli++)  
   {   

    //    Info<<qin_.internalField()[celli]<<endl; 
   } 
   
   forAll(mesh_.cells(),celli)  //for(label celli = 0; celli<mesh_.nCells(); celli++)  
   {

    //    Info<<qr_.internalField()[celli]<<endl;
   } 
}
void Foam::radiation::oqfvDOM::updateG()
{
    G_ = dimensionedScalar(dimMass/pow3(dimTime), Zero);
    qr_ = dimensionedScalar(dimMass/pow3(dimTime), Zero);
    qem_ = dimensionedScalar(dimMass/pow3(dimTime), Zero);
    qin_ = dimensionedScalar(dimMass/pow3(dimTime), Zero);

    forAll(IRay_, rayI)
    {
        IRay_[rayI].addIntensity();
        G_ += IRay_[rayI].I()*IRay_[rayI].omega();
        qr_.boundaryFieldRef() += IRay_[rayI].qr().boundaryField();
        qem_.boundaryFieldRef() += IRay_[rayI].qem().boundaryField();
        qin_.boundaryFieldRef() += IRay_[rayI].qin().boundaryField();
    }
   forAll(mesh_.cells(),celli)  //for(label celli = 0; celli<mesh_.nCells(); celli++)  
   {

    //    Info<<qr_.internalField()[celli]<<endl;//aded
   } 
}


void Foam::radiation::oqfvDOM::setRayIdLambdaId
(
    const word& name,
    label& rayId,
    label& lambdaId
) const
{
    // Assuming name is in the form: CHARS_rayId_lambdaId
    const auto i1 = name.find('_');
    const auto i2 = name.find('_', i1+1);

    rayId    = readLabel(name.substr(i1+1, i2-i1-1));
    lambdaId = readLabel(name.substr(i2+1));
}

void Foam::radiation::oqfvDOM::delqr()        //aded
{
    Gi =  dimensionedScalar("zero",dimMass/(dimLength*pow3(dimTime)), 0.0);  //- need to initialise becasue it take previous iteration value
    delqr_ =  dimensionedScalar("zero",dimMass/(dimLength*pow3(dimTime)), 0.0); //- need to initialise becasue it take previous iteration value
    
    for (label i=0; i < nLambda_; i++)
    {
        forAll(IRay_, rayI)
        {
             Gi += IRay_[rayI].ILambda(i)*IRay_[rayI].omega()*aLambda_[i];
        }
        
        delqr_ += ((aLambda_[i]*4*blackBody_.bLambda(i)) - Gi); //- pi(3.14) will not be included here because ....
    }
}


const Foam::solarCalculator& Foam::radiation::oqfvDOM::solarCalc() const
{
    return solarCalculator_();
}

// Foam::label Foam::radiation::oqfvDOM::nBands() const          //aded
// {

// 	Foam::Info<<" Value of nBands in abso "<<nLambda_<<"\n";  //info here means cout.
//     return nLambda_;//.absorptionEmission_->nBands();
// }
// ************************************************************************* //
