// Author: Soroosh Sanatkhani
// University of Pittsburgh
// Created: January 1, 2017
// Last Modified: August 24, 2021

/* Classes used in DistributedVelocityInletFvPatchVectorField.C is defined in DistributedVelocityInletFvPatchVectorField.H */

//////////////////////////////////////////////////////////////////////
// * * * * * DistributedVelocityInletFvPatchVectorField.H * * * * *//
//////////////////////////////////////////////////////////////////////
#ifndef DistributedVelocityInletFvPatchVectorField_H
#define DistributedVelocityInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{
class DistributedVelocityInletFvPatchVectorField:
    public fixedValueFvPatchVectorField
{
public:

   //- Runtime type information
   TypeName("DistributedVelocityInlet");
   // Constructors
        //- Construct from patch and internal field
        DistributedVelocityInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        DistributedVelocityInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        DistributedVelocityInletFvPatchVectorField
        (
            const DistributedVelocityInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        DistributedVelocityInletFvPatchVectorField
        (
            const DistributedVelocityInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new DistributedVelocityInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        DistributedVelocityInletFvPatchVectorField
        (
            const DistributedVelocityInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new DistributedVelocityInletFvPatchVectorField(*this, iF)
            );
        }
    // Member functions
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};
} // End namespace Foam
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
