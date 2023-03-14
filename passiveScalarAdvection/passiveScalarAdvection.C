// Author: Soroosh Sanatkhani
// University of Pittsburgh
// Created: January 1, 2017
// Last Modified: August 24, 2021

/* Please cite our paper if you find this repository useful:
Sanatkhani, S., Nedios, S., Menon, P. G., Saba, S. F., Jain, S. K., Federspiel, W. J., & Shroff, S. G. (2023).
Subject-specific factors affecting particle residence time distribution of left atrial appendage in atrial fibrillation:
A computational model-based study. Front Cardiovasc Med, 10(1070498), 1-13. https://doi.org/10.3389/fcvm.2023.1070498 */

//////////////////////////////////////////////////////////////////////
// * * * * * * *  * * * passiveScalarAdvection.C * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    simpleControl simple(mesh);
    #include "createFields.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nCalculating scalar transport\n" << endl;
    #include "CourantNo.H"
    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(s)
              + fvm::div(phi, s)
             ==
                fvOptions(s)
            );
            TEqn.relax();
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(s);
        }
        runTime.write();
    }
    Info<< "End\n" << endl;
    return 0;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
