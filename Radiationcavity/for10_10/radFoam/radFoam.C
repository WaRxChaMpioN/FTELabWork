/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    radFoam

Description
    Radiation solver only

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "radiationModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "scatterModel.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"  //it will be created formt he Control Dict file
    #include "createMesh.H"
    #include "createFvOptions.H"
  //  #include "initContinuityErrs.H"

	Info<< "Reading Temperature field\n" << endl;
	volScalarField T
	(
		IOobject
		(
			"T",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
	mesh
	);
	
	/*Info<< "Reading Absorption Coefficient field\n" << endl;
	volScalarField a
	(
		IOobject
		(
			"a",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		),
	mesh
	);*/

	autoPtr<radiation::radiationModel> radiation
	(
		radiation::radiationModel::New(T)
	);

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())  //Simple is the object which is created in the simpleControl.H and loop() func is from createfvOption.H. if we have loop(run time) then run time is the object for  the file setRoot Case.
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;  // this line will print the time reading from the Control Dict file and run the loop till the end time.
        {
        	radiation->correct();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
