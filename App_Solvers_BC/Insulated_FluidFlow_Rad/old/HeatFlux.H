     if (runTime.writeTime())
    {
     runTime.write();
     Info << "Reading field Qr\n" << endl;
    volScalarField qr
    (
        IOobject
        (
            "qr",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("qr", dimMass/(dimLength*dimTime*dimTime*dimTime),0)
    );

    
    surfaceScalarField gradT
    (
        IOobject
        (
            "gradT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("gradT", dimTemperature/dimLength,scalar(0))
    );

  
        gradT=fvc::snGrad(T);
        
        surfaceScalarField heatFlux =k*gradT;

        
        const surfaceScalarField::Boundary& patchGradT =gradT.boundaryField();
        
        const surfaceScalarField::Boundary& patchHeatFlux =heatFlux.boundaryField();
            
           scalar nt=-1;    
        const volScalarField Nhqr=nt*qr; 
        
        //for the printing of the output file only
        const volScalarField::Boundary& NegpatchRadHeatFlux = Nhqr.boundaryField();   
       
       volScalarField NegRadHeatFlux
        (
            IOobject
            (
                "NegRadHeatFlux",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("NegRadHeatFlux", qr.dimensions(), 0.0)
        );
 
          
            forAll(NegRadHeatFlux.boundaryField(), patchi)
            {
                NegRadHeatFlux.boundaryFieldRef()[patchi] =NegpatchRadHeatFlux[patchi];
                    
            }
    
          NegRadHeatFlux.write();   
             
             
             
        const volScalarField::Boundary& patchRadHeatFlux =  qr.boundaryField(); 
        
        const surfaceScalarField::Boundary& magSf = mesh.magSf().boundaryField();

        Info<< "\nWall heat fluxes [W]" << endl;
        forAll(patchHeatFlux, patchi)
        {
           // if (typeid(mesh.boundary()[patchi]) == typeid(wallFvPatch))
           // {
                scalar convFlux = gSum(magSf[patchi]*patchHeatFlux[patchi]);
                scalar radFlux = -gSum(magSf[patchi]*patchRadHeatFlux[patchi]);

                Info<< mesh.boundary()[patchi].name() << endl
                    << "    convective: " << convFlux << endl
                    << "    radiative:  " << radFlux << endl
                    << "    total:      " << convFlux + radFlux << endl;
          //  }
        }
        Info<< endl;

        volScalarField wallHeatFlux
        (
            IOobject
            (
                "wallHeatFlux",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("wallHeatFlux", heatFlux.dimensions(), 0.0)
        );

        forAll(wallHeatFlux.boundaryField(), patchi)
        {
            wallHeatFlux.boundaryFieldRef()[patchi] = patchHeatFlux[patchi];
        }

        wallHeatFlux.write();

        // Write the total heat-flux including the radiative contribution
        // if available
        if (qr.headerOk())
        {
            volScalarField totalWallHeatFlux
            (
                IOobject
                (
                    "totalWallHeatFlux",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar
                (
                    "totalWallHeatFlux",
                    heatFlux.dimensions(),
                    0.0
                )
            );

            forAll(totalWallHeatFlux.boundaryField(), patchi)
            {
                totalWallHeatFlux.boundaryFieldRef()[patchi] =
                    patchHeatFlux[patchi] - patchRadHeatFlux[patchi];
            }

            totalWallHeatFlux.write();
        }
        
        //calculation of the conduction Nu
        
        volScalarField CondNu
	  (
	     IOobject
		(
		   "CondNu",
		   runTime.timeName(),
		   mesh
	  	 ),
	     mesh,
	     dimensionedScalar("CondNu", heatFlux.dimensions(), 0.0)
	  );
        
       
                
        forAll(CondNu.boundaryField(), patchi)
	{
	  
	  CondNu.boundaryFieldRef()[patchi]=patchHeatFlux[patchi]/(20.6082);
	}
	
        CondNu.write();
        
       //calculation of the Radiation Nu
       volScalarField RadNu
	  (
	     IOobject
		(
		   "RadNu",
		   runTime.timeName(),
		   mesh
	  	 ),
	     mesh,
	     dimensionedScalar("RadNu",  dimless, 0.0)
	  );
        
        
        
        forAll(RadNu.boundaryField(), patchi)
	{
	  RadNu.boundaryFieldRef()[patchi]=-patchRadHeatFlux[patchi]/(20.6082);
	}
	
        RadNu.write();
   
       //calculation of the total Nu
        volScalarField TotalNu
	  (
	     IOobject
		(
		   "TotalNu",
		   runTime.timeName(),
		   mesh
	  	 ),
	     mesh,
	     dimensionedScalar("TotalNu",  dimless, 0.0)
	  );
        
       forAll(TotalNu.boundaryField(), patchi)
	{
	  TotalNu.boundaryFieldRef()[patchi]=((patchHeatFlux[patchi])/(20.6082))-
	                                       ((patchRadHeatFlux[patchi])/(20.6082));
	}
        TotalNu.write();
  /*      
       label patchI = mesh.boundaryMesh().findPatchID("leftWall"); 
       Info<< "The patch presently on "<<patchI<<endl;   
    
   
   forAll(mesh.boundaryMesh()[patchI].faceCentres(), faceI)
   {
    scalar x = mesh.boundaryMesh()[patchI].faceCentres()[faceI].x();
    scalar y = mesh.boundaryMesh()[patchI].faceCentres()[faceI].y();
    scalar z = mesh.boundaryMesh()[patchI].faceCentres()[faceI].z();
    Info<<faceI<<" "<<x<<" "<<y<<" "<<z<<" "<<endl;
   }     
*/
 }   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
