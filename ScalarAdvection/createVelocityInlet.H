// Author: Soroosh Sanatkhani
// University of Pittsburgh
// Created: January 1, 2017
// Last Modified: August 24, 2021

// Developed based on a sample createFields of OpenFOAM.

/* UniformVelocity variable array is defined to store the PV velocity inlets magnitude at each time step. */

//////////////////////////////////////////////////////////////////////
// * * * * * * * * * * createVelocityInlet.H * * * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
IOdictionary VelocityInletProperties
(
    IOobject
    (
        "VelocityInletProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const wordList titleNames(VelocityInletProperties.toc());

forAll(titleNames, item)
	{
    const word& titleName = titleNames[item];
    const dictionary& subDict = VelocityInletProperties.subDict(titleName);
	}

Info<< "Reading scalar list UniformVelocity \n" << endl;
scalarIOList UniformVelocity
(
    IOobject
    (
        "UniformVelocity",
        runTime.timeName(),
        mesh
    ),
    1
);
/*Initialising the storage array*/
for (int i = 0; i <1; i++)
{
  UniformVelocity[i] = 0;
}
/* Initialise the VelocityInlet Properties*/
initialiseVel(VelocityInletProperties); 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
