// Author: Soroosh Sanatkhani
// University of Pittsburgh
// Created: January 1, 2017
// Last Modified: August 24, 2021

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
