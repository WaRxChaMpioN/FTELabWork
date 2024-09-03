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

#include "fvDOM.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "constants.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(fvDOM, 0);
        addToRadiationRunTimeSelectionTables(fvDOM);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::fvDOM::initialise()
{
    // 3D
    if (mesh_.nSolutionD() == 3)
    {
        // Info<<mesh_.nSolutionD()<<endl;                                             //aded
        nRay_ = 4*nPhi_*nTheta_; //Total number of rays (1 per direction)
        IRay_.setSize(nRay_);  //Resize the vector
        scalar deltaPhi = pi/(2.0*nPhi_);  //dphi
        scalar deltaTheta = pi/nTheta_;  //dtheta
        label i = 0;  //label means int
        for (label n = 1; n <= nTheta_; n++)
        {
            for (label m = 1; m <= 4*nPhi_; m++)
            {
                scalar thetai = (2.0*n - 1.0)*deltaTheta/2.0;
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set  //set is the function which store value in the array in IRay_
                (
                    i,
                    new radiativeIntensityRay
                    (
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        absorptionEmission_,
                        scatter_,
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
        Info<<mesh_.nSolutionD()<<endl;             //aded
        scalar thetai = piByTwo;
        scalar deltaTheta = pi;
        nRay_ = 4*nPhi_;
        IRay_.setSize(nRay_);
        scalar deltaPhi = pi/(2.0*nPhi_);
        label i = 0;
        for (label m = 1; m <= 4*nPhi_; m++)
        {
            scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
            IRay_.set
            (
                i,
                new radiativeIntensityRay
                (
                    *this,
                    mesh_,
                    phii,
                    thetai,
                    deltaPhi,
                    deltaTheta,
                    nLambda_,
                    absorptionEmission_,
                    scatter_,
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
        Info<<mesh_.nSolutionD()<<endl;         //aded
        scalar thetai = piByTwo;
        scalar deltaTheta = pi;
        nRay_ = 2;
        IRay_.setSize(nRay_);
        scalar deltaPhi = pi;
        label i = 0;
        for (label m = 1; m <= 2; m++)
        {
            scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
            IRay_.set
            (
                i,
                new radiativeIntensityRay
                (
                    *this,
                    mesh_,
                    phii,
                    thetai,
                    deltaPhi,
                    deltaTheta,
                    nLambda_,
                    absorptionEmission_,
                    scatter_,
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
                    IOobject::AUTO_WRITE
                ),
                a_
            )
        );
    }

    Info<< "fvDOM : Allocated " << IRay_.size()
        << " rays with average orientation:" << nl;

    forAll(IRay_, rayId)
    {
        if (omegaMax_ <  IRay_[rayId].omega())
        {
            omegaMax_ = IRay_[rayId].omega();
        }
        Info<< '\t' << IRay_[rayId].I().name() << " : " << "dAve : "
            << '\t' << IRay_[rayId].dAve() << nl;
    }

    Info<< endl;

    if (this->found("useSolarLoad"))
    {
        this->lookup("useSolarLoad") >> useSolarLoad_;
    }

    if (useSolarLoad_)
    {
        /*const dictionary& solarDict = this->subDict("solarLoarCoeffs");
        solarLoad_.reset
        (
            new solarLoad(solarDict, T_, externalRadHeatFieldName_)
        );

        if (solarLoad_->nBands() > 1)
        {
            FatalErrorInFunction
                << "Requested solar radiation with fvDOM. Using "
                << "more than one band for the solar load is not allowed"
                << abort(FatalError);
        }

        Info<< "Creating Solar Load Model " << nl;
        */
    }

	//	const dictionary& nWeights
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::fvDOM::fvDOM(const volScalarField& T)
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
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
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
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qr", dimMass/pow3(dimTime), 0.0)
    ),
    qem_
    (
        IOobject
        (
            "qem",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qem", dimMass/pow3(dimTime), 0.0)
    ),
    qin_
    (
        IOobject
        (
            "qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qin", dimMass/pow3(dimTime), 0.0)
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
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    nTheta_(readLabel(coeffs_.lookup("nTheta"))),
    nPhi_(readLabel(coeffs_.lookup("nPhi"))),
    nRay_(0),
    nLambda_(readLabel(coeffs_.lookup("nWeights"))), //Uncomment it when not using the LBL case   //aded
    //nLambda_(absorptionEmission_->nBands()), // Uncomment it when using the LBL case   //removed
    aLambda_(nLambda_),
    blackBody_(nLambda_, T),
    IRay_(0),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    omegaMax_(0),
    useSolarLoad_(false),
    solarLoad_(),
    meshOrientation_
    (
        coeffs_.lookupOrDefault<vector>("meshOrientation", Zero)
    )
{
    initialise();
}


Foam::radiation::fvDOM::fvDOM
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
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
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
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qr", dimMass/pow3(dimTime), 0.0)
    ),
    qem_
    (
        IOobject
        (
            "qem",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qem", dimMass/pow3(dimTime), 0.0)
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
    qin_
    (
        IOobject
        (
            "qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qin", dimMass/pow3(dimTime), 0.0)
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
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    nTheta_(readLabel(coeffs_.lookup("nTheta"))),
    nPhi_(readLabel(coeffs_.lookup("nPhi"))),
    nRay_(0),
    nLambda_(readLabel(coeffs_.lookup("nWeights"))), //Uncomment it when not using the LBL case  //aded
    //nLambda_(readLabel(coeffs_.lookup("nBand"))), // Uncomment it when using the LBL case   //commented
    aLambda_(nLambda_),
    blackBody_(nLambda_, T),
    IRay_(0),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    omegaMax_(0),
    useSolarLoad_(false),
    solarLoad_(),
    meshOrientation_
    (
        coeffs_.lookupOrDefault<vector>("meshOrientation", Zero)
    )
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::fvDOM::~fvDOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::fvDOM::read()
{
    if (radiationModel::read())
    {
        // Only reading solution parameters - not changing ray geometry
        coeffs_.readIfPresent("convergence", convergence_);
        coeffs_.readIfPresent("maxIter", maxIter_);

        return true;
    }
    else
    {
        return false;
    }
}

void Foam::radiation::fvDOM::calculate()
{
    absorptionEmission_->correct(a_, aLambda_);

    updateBlackBodyEmission();

    if (useSolarLoad_)
    {
        solarLoad_->calculate();
    }
    
    // Set rays convergence false
    List<bool> rayIdConv(nRay_, false);

    scalar maxResidual = 0.0;
    label radIter = 0;
    do
    {
        Info<< "Radiation solver iter: " << radIter << endl; 
        //Do not delete    
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

                if (maxBandResidual < convergence_)
                {
                    rayIdConv[rayI] = true;
                }
            }
        }
    } while (maxResidual > convergence_ && radIter < maxIter_);

    updateG();
    delqr();            //aded
    qr_Domain();           //aded 
    qr_Boundary();      //aded
    energyCheck();      //aded
}


Foam::tmp<Foam::volScalarField> Foam::radiation::fvDOM::Rp() const// I wanted to print this function that is problem because many paper report del.qr which is Ru and Rp()
{
  

  List<scalar> aj (coeffs_.lookup("weight"));   //aded

    // Construct using contribution from first frequency band
    tmp<volScalarField> tRp
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
			    false
                
            ),
            (
                4
               *physicoChemical::sigma
               *(aLambda_[0] - absorptionEmission_->aDisp(0)())     //aded
               *blackBody_.deltaLambdaT(T_, absorptionEmission_->bands(0)) //aded
									
            )
						
        )
    );
	
    volScalarField& Rp=tRp.ref();

    // Add contributions over remaining frequency bands
    for (label j=1; j < nLambda_; j++)
    {
        Rp +=
        (
            4
           *physicoChemical::sigma
           *(aLambda_[j] - absorptionEmission_->aDisp(j)())
           *blackBody_.deltaLambdaT(T_, absorptionEmission_->bands(j))*aj[j]
        );

	
    }
		


    return tRp;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiation::fvDOM::Ru() const   //This function too Ru-Rp is our source term del.qr fro energy equation, thath we find in radiation model Ru was earlier in this 
{
                            // added whole box
     List<scalar> aj (coeffs_.lookup("weight"));

    tmp<DimensionedField<scalar, volMesh>> tRu
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "Ru",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
				false
                
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0)
        )
    );

    DimensionedField<scalar, volMesh>& Ru=tRu.ref();

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
        Ru += (aLambda_[j]() - absorptionEmission_->aDisp(j)()())*Gj*aj[j]
             - absorptionEmission_->ECont(j)()();	
    }

    return tRu;
}


void Foam::radiation::fvDOM::updateBlackBodyEmission()
{
    for (label j=0; j < nLambda_; j++)
    {
        blackBody_.correct(j, absorptionEmission_->bands(j));
    }
}

void Foam::radiation::fvDOM::qr_Domain()            //aded
{
   qr_Domain_ = 0.0; //- need to initialise becasue it take previous iteration value
   
   forAll(mesh_.cells(),celli)  //for(label celli = 0; celli<mesh_.nCells(); celli++)  
   {
       qr_Domain_ += delqr_[celli]*mesh_.V()[celli];
   } 
}

void Foam::radiation::fvDOM::qr_Boundary()          //aded
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

void Foam::radiation::fvDOM::energyCheck()              //aded
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

void Foam::radiation::fvDOM::qr_Cell()                  //aded
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



void Foam::radiation::fvDOM::updateG()
{
    G_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
    qr_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
    qem_ = dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);
    qin_ = dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);

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



void Foam::radiation::fvDOM::setRayIdLambdaId
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

void Foam::radiation::fvDOM::delqr()        //aded
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



/*void Foam::radiation::fvDOM::delqr() // this is member function I defined to write del.qr 
{
                List<scalar> aj (coeffs_.lookup("weight"));
                const volScalarField T3(pow3(T_));


				qRu = dimensionedScalar("zero",dimMass/(dimLength*pow3(dimTime)), 0.0);

				qRP = dimensionedScalar("zero",dimMass/(dimLength*pow3(dimTime)), 0.0);

				Gi = dimensionedScalar("zero",dimMass/(dimLength*pow3(dimTime)), 0.0);


	forAll(mesh_.V(),celli)	
   	{
		

        for (label j=0; j < nLambda_; j++)
        {
				
		for (label rayI=0; rayI < nRay_; rayI++)
        {
            Gi[celli] += IRay_[rayI].ILambda(j)[celli]*IRay_[rayI].omega()*aLambda_[j][celli]*aj[j]; //calculates G and aj is weight function
        }

         qRu[celli] = (aLambda_[j][celli])*Gi[celli];  // this calcultes irradiation 
				
				qRP[celli] += (aLambda_[j][celli])*4*5.67037441918442945397*T3[celli]*T_[celli]*(1e-8)*aj[j];// this function is rough code I will fine tune it. and aj is weight function

                // this calculates del.qr
		}
			
		}

delq = (qRP-Gi);

}*/
/**************************************************************************************/

Foam::label Foam::radiation::fvDOM::nBands() const          //aded
{

	Foam::Info<<" Value of nBands in abso "<<nLambda_<<"\n";  //info here means cout.
    return nLambda_;//.absorptionEmission_->nBands();
}










// ************************************************************************* //
