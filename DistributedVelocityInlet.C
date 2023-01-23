// Author: Soroosh Sanatkhani
// University of Pittsburgh
// Created: January 1, 2017
// Last Modified: August 24, 2021

//////////////////////////////////////////////////////////////////////
// * * * * * * * * * DistributedVelocityInlet.C * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
double TotalArea;
std::string TotalFlowRateFile;

void initialiseVel(const dictionary& VelocityInletProperties)
{
 /* Initialising */

	const wordList titleNames(VelocityInletProperties.toc());

	forAll(titleNames, item)
	{
		const word& titleName = titleNames[item];

		const dictionary& subDict = VelocityInletProperties.subDict(titleName);

		TotalFlowRateFile = string(subDict.lookup("TotalFlowRateFile"));
   }
}

void Calculate_TotalArea(fvMesh & mesh)
{
	label pv1 = mesh.boundaryMesh().findPatchID("pv1");
	label pv2 = mesh.boundaryMesh().findPatchID("pv2");
	label pv3 = mesh.boundaryMesh().findPatchID("pv3");
	label pv4 = mesh.boundaryMesh().findPatchID("pv4");
	TotalArea =		gSum(mesh.magSf().boundaryField()[pv1])
				+	gSum(mesh.magSf().boundaryField()[pv2])
				+	gSum(mesh.magSf().boundaryField()[pv3])
				+	gSum(mesh.magSf().boundaryField()[pv4]);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
