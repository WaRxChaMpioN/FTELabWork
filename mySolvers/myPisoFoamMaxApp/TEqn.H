
forAll(Keff, celli)
{
    // Correctly access the face area (magnitude of the surface vector field)
    scalar Sf = (mesh.magSf()[celli]);
    
    // Correctly extract the time step size as a scalar
    scalar dt = runTime.deltaT().value();
    
    // Correctly access the cell volume
    scalar Vol = mesh.V()[celli];
    
    // Calculate the magnitude of the average velocity for the cell
    scalar averageVelocity = (U[celli].x());
    
    // Calculate the volume fraction v
    scalar v = (averageVelocity * Sf * dt) / Vol;
   	
	// thermal conductivity
	scalar Kf =Kfluid.value();
	scalar Ks = Ksolid.value();
	
	// numerator
	scalar num = Kf + 2*Ks - 2*v*(Ks - Kf);
	
	//Deniminator
	scalar den = Kf + 2*Ks + v*(Ks - Kf);

    // Calculate the effective thermal conductivity
	Keff[celli] = (num/den)*Ks;

}

volScalarField DT = Keff/(rho*CP);

// Temperature equation
fvScalarMatrix TEqn
(
   // fvm::ddt(T) +
     fvm::div(phi, T)
    - fvm::laplacian(Keff/(rho*CP), T)  // Use Keff in the laplacian term
    
    // Include porosity effects in the temperature equation
);

// Solve the temperature equation
solve(TEqn);

// Correct for any additional effects from fvOptions
//fvOptions.correct(T);
