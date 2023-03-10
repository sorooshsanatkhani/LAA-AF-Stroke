// Author: Soroosh Sanatkhani
// University of Pittsburgh
// Created: January 1, 2017
// Last Modified: August 24, 2021

// Developed based on a sample createFields of OpenFOAM.

//////////////////////////////////////////////////////////////////////
// * * * * * * * * * * * * * * * createFields.H * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

singlePhaseTransportModel fluid(U, phi);


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
