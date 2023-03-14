// Author: Soroosh Sanatkhani
// University of Pittsburgh
// Created: January 1, 2017
// Last Modified: August 24, 2021

/* Please cite our paper if you find this repository useful:
Sanatkhani, S., Nedios, S., Menon, P. G., Saba, S. F., Jain, S. K., Federspiel, W. J., & Shroff, S. G. (2023).
Subject-specific factors affecting particle residence time distribution of left atrial appendage in atrial fibrillation:
A computational model-based study. Front Cardiovasc Med, 10(1070498), 1-13. https://doi.org/10.3389/fcvm.2023.1070498 */

//////////////////////////////////////////////////////////////////////
// * * * * * * * * * * * * * * * Quemada.H * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
#ifndef Quemada_H
#define Quemada_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*------------------------------------------------------------------*\
                           Class Quemada Declaration
\*------------------------------------------------------------------*/
class Quemada:
    public viscosityModel
{
    // Private Data

        dictionary QuemadaCoeffs_;

        dimensionedScalar nuP_;
        dimensionedScalar nuMax_;
        dimensionedScalar H_;

        volScalarField nu_;

    // Private Member Functions
        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;

public:
    //- Runtime type information
    TypeName("Quemada");

    // Constructors
        //- Construct from components
        Quemada
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );
    //- Destructor
    virtual ~Quemada()
    {}

    // Member Functions
        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }
        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace viscosityModels
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
