/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fluidPhaseModel.H"
#include "twoPhaseSystem.H"
#include "diameterModel.H"
#include "fvMatrix.H"
#include "dragModel.H"
#include "heatTransferModel.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        fluidPhaseModel,
        dictionary,
        fluid
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidPhaseModel::fluidPhaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    phaseModel(fluid, phaseProperties, phaseName),
    p_
    (
        IOobject
        (
            IOobject::groupName("p", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->thermoPtr_->p(),
        this->thermoPtr_->p().boundaryField().types()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidPhaseModel::~fluidPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluidPhaseModel::advect
(
    const dimensionedScalar& deltaT,
    const dimensionedVector& g,
    const volVectorField& Ui,
    const volScalarField& pi,
    const bool oldTime
)
{
    volScalarField& alpha = *this;
    if (oldTime)
    {
        if (!otherPhase().granular())
        {
            alpha = alpha.oldTime() - deltaT*(Ui & gradAlpha_);
            alpha.correctBoundaryConditions();
        }

        alphaRho_ = alphaRho_.oldTime() - deltaT*fvc::div(massFlux_);
        alphaRho_.correctBoundaryConditions();

        alphaRhoU_ =
            alphaRhoU_.oldTime()
          + deltaT
           *(
              - fvc::div(momentumFlux_)
              + p_*gradAlpha_
              + alpha*rho_*g
            );
        alphaRhoU_.correctBoundaryConditions();

        alphaRhoE_ =
            alphaRhoE_.oldTime()
          + deltaT
           *(
              - fvc::div(energyFlux_)
              - pi*fvc::div(otherPhase().massFlux())/otherPhase().rho()
              + (alpha*rho_*g & U_)
            );
        alphaRhoE_.correctBoundaryConditions();
    }
    else
    {
        if (!otherPhase().granular())
        {
            alpha -= deltaT*(Ui & gradAlpha_);
            alpha.correctBoundaryConditions();
        }

        alphaRho_ -= deltaT*fvc::div(massFlux_);
        alphaRho_.correctBoundaryConditions();

        alphaRhoU_ +=
            deltaT
           *(
              - fvc::div(momentumFlux_)
              + p_*gradAlpha_
              + alpha*rho_*g
            );
        alphaRhoU_.correctBoundaryConditions();

        alphaRhoE_ +=
            deltaT
           *(
              - fvc::div(energyFlux_)
              - pi*fvc::div(otherPhase().massFlux())/otherPhase().rho()
              + (alpha*rho_*g & U_)
            );
        alphaRhoE_.correctBoundaryConditions();
    }
}


void Foam::fluidPhaseModel::solveSources
(
    const dimensionedScalar& deltaT,
    const volScalarField& Mdot,
    const volVectorField& Ftot,
    const volScalarField& qConv,
    const volScalarField& gamma,
    const volVectorField& Ui,
    const volScalarField& pi
)
{
    alphaRho_ += deltaT*Mdot;
    alphaRhoU_ += deltaT*Ftot;
    alphaRhoE_ += deltaT*(gamma + qConv + (Ftot & otherPhase().U()));
}


void Foam::fluidPhaseModel::updateFluxes()
{
    volScalarField H(IOobject::groupName("H", name()), E_ + p_/rho_);
    if (otherPhase().granular())
    {
        this->fluxFunction_->updateFluxes
        (
            alphaf_,
            massFlux_,
            momentumFlux_,
            energyFlux_,
            *this,
            rho_,
            U_,
            H,
            p_,
            c(),
            fluid_.U(),
            fluid_.p()
        );
    }
    else
    {
        this->fluxFunction_->updateFluxes
        (
            massFlux_,
            momentumFlux_,
            energyFlux_,
            alphaf_,
            rho_,
            U_,
            H,
            p_,
            c(),
            fluid_.U(),
            fluid_.p()
        );
    }
    gradAlpha_ = fvc::surfaceIntegrate(fluid_.mesh().Sf()*alphaf_);
}


void Foam::fluidPhaseModel::decode()
{
    (*this).max(residualAlpha_);
    (*this).correctBoundaryConditions();

    rho_ = alphaRho_/(*this);
    rho_.correctBoundaryConditions();

    U_ = alphaRhoU_/alphaRho_;
    phiPtr_() = fvc::flux(U_);

    E_ = alphaRhoE_/alphaRho_;
    he_ = E_ - 0.5*magSqr(U_);

    thermoPtr_->correctP();
    p_ = thermoPtr_->p();
    p_.correctBoundaryConditions();
}


void Foam::fluidPhaseModel::encode()
{
    (*this).max(residualAlpha_);
    (*this).correctBoundaryConditions();

    alphaRho_ = (*this)*rho_;
    alphaRho_.correctBoundaryConditions();

    alphaRhoU_ = alphaRho_*U_;
    alphaRhoU_.correctBoundaryConditions();

    E_ = he_ + 0.5*magSqr(U_);
    alphaRhoE_ = alphaRho_*E_;
    alphaRhoE_.correctBoundaryConditions();
}


void Foam::fluidPhaseModel::correctThermo()
{
    E_ = he_ + 0.5*magSqr(U_);

    thermoPtr_->correctP();
    p_ = thermoPtr_->p();
    p_.correctBoundaryConditions();
}

void Foam::fluidPhaseModel::store()
{
    (*this).storeOldTime();
    alphaRho_.storeOldTime();
    alphaRhoU_.storeOldTime();
    alphaRhoE_.storeOldTime();
}

// ************************************************************************* //
