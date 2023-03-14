// Author: Soroosh Sanatkhani
// University of Pittsburgh
// Created: January 1, 2017
// Last Modified: August 24, 2021

/* Please cite our paper if you find this repository useful:
Sanatkhani, S., Nedios, S., Menon, P. G., Saba, S. F., Jain, S. K., Federspiel, W. J., & Shroff, S. G. (2023).
Subject-specific factors affecting particle residence time distribution of left atrial appendage in atrial fibrillation:
A computational model-based study. Front Cardiovasc Med, 10(1070498), 1-13. https://doi.org/10.3389/fcvm.2023.1070498 */

// Developed based on icoFoam of OpenFOAM.

//////////////////////////////////////////////////////////////////////
// * * * * * * * * * * ScalarAdvection.C * * * * * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
#include "fvCFD.H"
#include "pisoControl.H"
#include "scalarIOList.H"
#include "DistributedVelocityInlet.C"
#include "DistributedVelocityInletFvPatchVectorField.H"

///////// For CSV Reader //////////////
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createVelocityInlet.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

//////////////////////////////////////////////////////////////////////
    std::vector<std::vector<std::vector<double>>> FlowRate;
    for (int i = 0; i<1; i++)
    {
    	std::ifstream in(TotalFlowRateFile);
    	std::vector<std::vector<double>> FlowRate2;
    	if (in)
    	{
    		std::string line;
    		while (getline(in, line))
    		{
    			std::stringstream sep(line);
    			std::string flowrate;
    			FlowRate2.push_back(std::vector<double>());
    			while (getline(sep, flowrate, ','))
    			{
    				FlowRate2.back().push_back(stod(flowrate));
    			}
    		}
    	}
    	FlowRate.push_back(FlowRate2);
    }

    Calculate_TotalArea(mesh);
  //////////////////////////////////////////////////////////////////////

    Info<< "\nStarting time loop\n" << endl;
	
	int d = int(round(1000*runTime.value())) % FlowRate[0].size();
	
    UniformVelocity[0] = (FlowRate[0][d][1]/60000)/TotalArea;

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

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        solve(fvm::ddt(s) + fvm::div(phi,s));
		
		#include "continuityErrs.H"

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
