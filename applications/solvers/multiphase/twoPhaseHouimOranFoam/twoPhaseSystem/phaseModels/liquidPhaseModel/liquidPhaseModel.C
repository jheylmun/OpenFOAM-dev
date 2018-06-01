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

#include "liquidPhaseModel.H"
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
        liquidPhaseModel,
        dictionary,
        liquid
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquidPhaseModel::liquidPhaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    phaseModel(fluid, phaseProperties, phaseName),
    p_(otherPhase().p())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liquidPhaseModel::~liquidPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::liquidPhaseModel::advect
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
        alpha = alpha.oldTime() - deltaT*(Ui & gradAlpha_);
        alpha.correctBoundaryConditions();

        this->alphaRho_ =
            this->alphaRho_.oldTime()
          - deltaT*fvc::div(this->massFlux_);
        this->alphaRho_.correctBoundaryConditions();

        this->alphaRhoU_ =
            this->alphaRhoU_.oldTime();
          - deltaT
           *(
                fvc::div(this->momentumFlux_)
              - pi*this->gradAlpha_
              + alpha*this->rho()*g
            );
        this->alphaRhoU_.correctBoundaryConditions();

        this->alphaRhoE_ =
            this->alphaRhoE_.oldTime()
          - deltaT
           *(
                fvc::div(energyFlux_)
              - pi*(Ui & this->gradAlpha_)
              + (alpha*this->rho()*g & this->U_)
            );
        this->alphaRhoE_.correctBoundaryConditions();
    }
    else
    {
        alpha -= deltaT*(Ui & gradAlpha_);
        alpha.correctBoundaryConditions();

        this->alphaRho_ -= deltaT*fvc::div(this->massFlux_);
        this->alphaRho_.correctBoundaryConditions();

        this->alphaRhoU_ -=
            deltaT
           *(
                fvc::div(this->momentumFlux_)
              - pi*this->gradAlpha_
              + alpha*this->rho()*g
            );
        this->alphaRhoU_.correctBoundaryConditions();

        this->alphaRhoE_ -=
            deltaT
           *(
                fvc::div(energyFlux_)
              - pi*(Ui & this->gradAlpha_)
              + (alpha*this->rho()*g & this->U_)
            );
        this->alphaRhoE_.correctBoundaryConditions();
    }
}


void Foam::liquidPhaseModel::solveSources
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
    alphaRhoE_ += deltaT*(gamma - qConv);
}

void Foam::liquidPhaseModel::updateFluxes()
{
    // calculate fluxes with
    volScalarField H(IOobject::groupName("H", name()), E_ + p_/rho_);
    this->fluxFunction_->updateFluxes
    (
        gradAlpha_,
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


void Foam::liquidPhaseModel::decode()
{
    (*this).max(this->residualAlpha_);
    (*this).correctBoundaryConditions();

    this->rho_ = alphaRho_/(*this);
    this->rho_.correctBoundaryConditions();

    this->U_ = this->alphaRhoU_/this->alphaRho_;
    this->U_.correctBoundaryConditions();
    this->phiPtr_() = fvc::flux(this->U_);

    this->E_ = this->alphaRhoE_/this->alphaRho_;
    this->he_ -= 0.5*magSqr(this->U_);
}


void Foam::liquidPhaseModel::encode()
{
    (*this).max(this->residualAlpha_);
    (*this).correctBoundaryConditions();

    this->alphaRho_ = (*this)*this->rho_;
    this->alphaRho_.correctBoundaryConditions();

    this->alphaRhoU_ = this->alphaRho_*this->U_;
    this->alphaRhoU_.correctBoundaryConditions();

    this->alphaRhoE_ = this->alphaRho_*this->E_;
    this->alphaRhoE_.correctBoundaryConditions();
}


void Foam::liquidPhaseModel::correctThermo()
{
    this->E_ = this->he_ + 0.5*magSqr(U_);
}

void Foam::liquidPhaseModel::store()
{
    (*this).storeOldTime();
    alphaRho_.storeOldTime();
    alphaRhoU_.storeOldTime();
    alphaRhoE_.storeOldTime();
}

// ************************************************************************* //
