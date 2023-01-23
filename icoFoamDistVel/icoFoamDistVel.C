// Author: Soroosh Sanatkhani
// University of Pittsburgh
// Created: January 1, 2017
// Last Modified: August 24, 2021

//////////////////////////////////////////////////////////////////////
// * * * * * * * * * * * * * * icoFoamDistVel.C * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
#include "fvCFD.H"
#include "pisoControl.H"
#include "scalarIOList.H"
#include "DistributedVelocityInlet.C"
#include "DistributedVelocityInletFvPatchVectorField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    pisoControl piso(mesh);
    #include "createVelocityInlet.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	Calculate_TotalArea(mesh);
	UniformVelocity[0] = (TotalFlowRate/60000)/TotalArea;
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
                pEqn.solve();
                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }
            #include "continuityErrs.H"
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }
        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "End\n" << endl;
    return 0;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
