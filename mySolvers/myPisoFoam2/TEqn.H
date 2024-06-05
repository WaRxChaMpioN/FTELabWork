
forAll(Keff, celli)
{
    // Correctly access the face area (magnitude of the surface vector field)
    scalar Sf = mag(mesh.Sf()[celli]);
    
    // Correctly extract the time step size as a scalar
    scalar dt = runTime.deltaT().value();
    
    // Correctly access the cell volume
    scalar Vol = mesh.V()[celli];
    
    // Calculate the magnitude of the average velocity for the cell
    scalar averageVelocity = mag(U[celli]);
    
    // Calculate the volume fraction v
    scalar v = (averageVelocity * Sf * dt) / Vol;
    
    // Calculate the effective thermal conductivity
	Keff[celli] = v * Kfluid.value() + (1 - v) * Ksolid.value();


}

volScalarField DT = Keff/(rho*CP);

// Temperature equation
fvScalarMatrix TEqn
(
    fvm::ddt(T) +
     fvm::div(phi, T)
    - fvm::laplacian(Keff/(rho*CP), T)  // Use Keff in the laplacian term
    
    // Include porosity effects in the temperature equation
);

// Solve the temperature equation
solve(TEqn);

// Correct for any additional effects from fvOptions
fvOptions.correct(U);
