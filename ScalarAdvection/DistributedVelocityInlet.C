// Author: Soroosh Sanatkhani
// University of Pittsburgh
// Created: January 1, 2017
// Last Modified: August 24, 2021

/* Please cite our paper if you find this repository useful:
Sanatkhani, S., Nedios, S., Menon, P. G., Saba, S. F., Jain, S. K., Federspiel, W. J., & Shroff, S. G. (2023).
Subject-specific factors affecting particle residence time distribution of left atrial appendage in atrial fibrillation:
A computational model-based study. Front Cardiovasc Med, 10(1070498), 1-13. https://doi.org/10.3389/fcvm.2023.1070498 */

/* The total PV inlets surface area needed for UniformVelocity calculation is computed in DistributedVelocityInlet.C */

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
