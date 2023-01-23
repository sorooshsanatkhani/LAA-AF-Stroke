// Author: Soroosh Sanatkhani
// University of Pittsburgh
// Created: January 1, 2017
// Last Modified: August 24, 2021

//////////////////////////////////////////////////////////////////////
// * * * * * * * * * * * * * * * Quemada.C * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
#include "Quemada.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
// * * * * * * * * * * Static Data Members * * * * * * * * * * * * //
namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(Quemada, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        Quemada,
        dictionary
    );
}
}
// * * * * * * * * Private Member Functions  * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::Quemada::calcNu() const
{
	return min(
	nuMax_,
	nuP_*pow(1.0-0.5*
	((exp(3.874 - 10.41*H_ + 13.8*H_*H_ - 6.738*H_*H_*H_)
	+exp(1.3435 - 2.803*H_ + 2.711*H_*H_ - 0.6479*H_*H_*H_)
	*sqrt((dimensionedScalar(dimTime, 1.0)*strainRate())/exp(-6.1508 + 27.923*H_ - 25.6*H_*H_ + 3.697*H_*H_*H_)))
	/(1.0+sqrt((dimensionedScalar(dimTime, 1.0)*strainRate())/exp(-6.1508 + 27.923*H_ - 25.6*H_*H_ + 3.697*H_*H_*H_))))
	*H_
	,-2));
}
// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * //
Foam::viscosityModels::Quemada::Quemada
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    QuemadaCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    nuP_("nuP", dimViscosity, QuemadaCoeffs_),
    nuMax_("nuMax", dimViscosity, QuemadaCoeffs_),
    H_("H", dimless, QuemadaCoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}
// * * *  * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::viscosityModels::Quemada::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    QuemadaCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    QuemadaCoeffs_.lookup("nuP") >> nuP_;
    QuemadaCoeffs_.lookup("nuMax") >> nuMax_;
    QuemadaCoeffs_.lookup("H") >> H_;

    return true;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
//////////////////////////////////////////////////////////////////////
