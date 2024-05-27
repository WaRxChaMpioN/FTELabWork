/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Application
    icoFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

    \heading Solver details
    The solver uses the PISO algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}}
          + \div \left( \vec{U} \vec{U} \right)
          - \div \left(\nu \grad \vec{U} \right)
          = - \grad p
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, laminar flow"
        " of Newtonian fluids."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }
        
/********************************************************/
//Temperature equation:
        fvScalarMatrix TEqn
        (
            fvm::ddt(T)
          + fvm::div(phi, T)
          - fvm::laplacian(DT, T)
        );
   solve(TEqn);
  
/********************************************************/
                
	runTime.write();
       runTime.printExecutionTime(Info);
	
    }
    
    /********************************************************/
    
    //calculation of bulk Mean Temperature
        

        int x_=20;
	int z_=350;
//	grad = Foam::fvc::snGrad(T).BoudarySurfield(j)
	//double Temp = 293;
	double Dia = 0.002;
	double dx = 0.5*Dia/x_ ;
	double Um = 0.011;
	//double q = 1200;
        List<scalar> T_mean(mesh.nCells()/x_);
        List<scalar> Nu(mesh.nCells()/x_);
      //  List<scalar> gradT = Foam::fvc::grad(T);
/*        
	for (int i = 0; i <(z_); i++)
	{
		double T_sum = 0;
		int j = i*z_;
		while (j < (x_*z_))
		{
		
				Info<<T[i]<<endl;
			T_sum += T[j];
			j=j+300;		
		}
		
		T_mean[i]= T_sum/x_;	
		Info<<T_mean[i]<<endl;
		
	}
	*/
	
	
		int j = 0;	
	
	for (int i = 0; i <(z_); i++)
	{
	//	int k = x_;
		double T_sum = 0;
		double gradX = 0;
		int iter =1;
		while (j < (i+1)*x_)
		{	//double temp = sqr(U[j] .x()) +sqr(U[j] .y()) +sqr(U[j] .z());
			//double Ures = Foam::sqrt(temp);
			//double temp = U[j];
			//Info<<j<<endl;			
			T_sum += (U[j] .z())*T[j]*(iter)*dx*dx;
			//Info<< j<< " "<< U[j] .z()<< " "<<T_sum<<endl;
			j++;
//			k--;
			iter++;		
		}
		
		
		
		T_mean[i] = 8*T_sum/(Um*Dia*Dia);
		
		
		//Info<<endl;
		int iterm = x_*(i+1)-1;
		int iterm1 = (x_)*(i+1)-2;
		
		T_mean[i] = 8*T_sum/(Um*Dia*Dia);
		
		gradX = 2*(T[iterm] - T[iterm1])/(Dia);
		//tmp<GeometricField<double, fvPatchField, surfaceMesh> > snGradT = fvc::snGrad(T);
		//vector gradient = snGradT()[j];
		//double gradX = gradient[0];
		Nu[i] = 0.35* ((gradX)/(-T[iterm] + T_mean[i]) );
		
		
		//H[i] = q/(T[20*i+19]-T_mean[i]);
		
		//Info << T_mean[i]<< endl;
		Info << Nu[i]<< endl;
	

	/*
	int iterm = x_*(i+1)-1;
	int iterm1 = (x_)*(i+1)-2;
	
	Info << T[iterm1]<< endl;
	Info << T[iter]<< endl;
	
	Nu = 2*x_*(T[iterm]-T[iterm1])/(T[iterm] - T_mean[i]);
		
		Info<<T_mean[i]<<" "<< Nu[i] << endl;
	
	
	*/
	
    

    }
 //   Info << endl << dx;
   // Info<<mesh.nCells();
    
/********************************************************/
    


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
