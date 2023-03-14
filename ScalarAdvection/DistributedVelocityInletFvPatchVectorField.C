// Author: Soroosh Sanatkhani
// University of Pittsburgh
// Created: January 1, 2017
// Last Modified: August 24, 2021

/* Please cite our paper if you find this repository useful:
Sanatkhani, S., Nedios, S., Menon, P. G., Saba, S. F., Jain, S. K., Federspiel, W. J., & Shroff, S. G. (2023).
Subject-specific factors affecting particle residence time distribution of left atrial appendage in atrial fibrillation:
A computational model-based study. Front Cardiovasc Med, 10(1070498), 1-13. https://doi.org/10.3389/fcvm.2023.1070498 */

/* Velocity vector field normal to each individual PV inlet is calculated in DistributedVelocityInletFvPatchVectorField.C */

//////////////////////////////////////////////////////////////////////
// * * * * * DistributedVelocityInletFvPatchVectorField.C * * * * *//
//////////////////////////////////////////////////////////////////////

#include "DistributedVelocityInletFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"
#include "vectorField.H"
#include "fvc.H"
#include "scalarIOList.H"

// * * * * ** * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::DistributedVelocityInletFvPatchVectorField::
DistributedVelocityInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
 (Abecasis et al., 2009)

Foam::DistributedVelocityInletFvPatchVectorField::
DistributedVelocityInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict)
{}

Foam::DistributedVelocityInletFvPatchVectorField::
DistributedVelocityInletFvPatchVectorField
(
    const DistributedVelocityInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper   
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{}

Foam::DistributedVelocityInletFvPatchVectorField::
DistributedVelocityInletFvPatchVectorField
(const DistributedVelocityInletFvPatchVectorField& ptf):
    fixedValueFvPatchField<vector>(ptf)
{}

Foam::DistributedVelocityInletFvPatchVectorField::
DistributedVelocityInletFvPatchVectorField
(
    const DistributedVelocityInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{}
// * * * * * * ** * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DistributedVelocityInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    /* Accessing the variables stored in mesh */
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const scalarIOList& UniformVelocity = mesh.lookupObject<scalarIOList>("UniformVelocity");

    const scalar U = UniformVelocity[0]; 
	operator==(-U*patch().nf());
	
    fixedValueFvPatchField<vector>::updateCoeffs();
}

void Foam::DistributedVelocityInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
}
namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       DistributedVelocityInletFvPatchVectorField
   );
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
