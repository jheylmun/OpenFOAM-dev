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
    p_(thermoPtr_->p())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liquidPhaseModel::~liquidPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::liquidPhaseModel::advect
(
    const label stepi,
    const scalarList& coeffs,
    const scalarList& Fcoeffs,
    const dimensionedScalar& deltaT,
    const dimensionedVector& g,
    const volVectorField& F,
    const volVectorField& Ui,
    const volScalarField& pi
)
{
    alphaRhos_[stepi] = alphaRho_;
    alphaRhoUs_[stepi] = alphaRhoU_;
    alphaRhoEs_[stepi] = alphaRhoE_;

    deltaAlphaRho_[stepi] = -fvc::div(massFlux_);
    deltaAlphaRhoU_[stepi] = -fvc::div(momentumFlux_) + alphaRho_*g + F;
    deltaAlphaRhoE_[stepi] = -fvc::div(energyFlux_) + ((*this)*rho_*g & U_);

    volScalarField alphaRho(alphaRho_*coeffs[stepi]);
    volVectorField alphaRhoU(alphaRhoU_*coeffs[stepi]);
    volScalarField alphaRhoE(alphaRhoE_*coeffs[stepi]);

    volScalarField deltaAlphaRho(deltaAlphaRho_[stepi]*Fcoeffs[stepi]);
    volVectorField deltaAlphaRhoU(deltaAlphaRhoU_[stepi]*Fcoeffs[stepi]);
    volScalarField deltaAlphaRhoE(deltaAlphaRhoE_[stepi]*Fcoeffs[stepi]);

    for (label i = 0; i < stepi; i++)
    {
        alphaRho += alphaRhos_[i]*coeffs[i];
        alphaRhoU += alphaRhoUs_[i]*coeffs[i];
        alphaRhoE += alphaRhoEs_[i]*coeffs[i];

        deltaAlphaRho += deltaAlphaRho_[i]*Fcoeffs[i];
        deltaAlphaRhoU += deltaAlphaRhoU_[i]*Fcoeffs[i];
        deltaAlphaRhoE += deltaAlphaRhoE_[i]*Fcoeffs[i];
    }

    alphaRho_ = alphaRho + deltaT*deltaAlphaRho;
    alphaRho_.correctBoundaryConditions();

    alphaRhoU_ = alphaRhoU + deltaT*deltaAlphaRhoU;
    alphaRhoU_.correctBoundaryConditions();

    alphaRhoE_ = alphaRhoE + deltaT*deltaAlphaRhoE;
    alphaRhoE_.correctBoundaryConditions();
}


void Foam::liquidPhaseModel::updateFluxes()
{
    volScalarField H(IOobject::groupName("H", name()), E_ + p_/rho_);
    this->fluxFunction_->updateFluxes
    (
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
    gradAlpha_ = fluxFunction_->gradAlpha();
}

void Foam::liquidPhaseModel::updateFluxes(const surfaceScalarField& alphaf)
{
    volScalarField H(IOobject::groupName("H", name()), E_ + p_/rho_);
    this->fluxFunction_->updateFluxes
    (
        massFlux_,
        momentumFlux_,
        energyFlux_,
        alphaf,
        rho_,
        U_,
        H,
        p_,
        c(),
        fluid_.U(),
        fluid_.p()
    );
    gradAlpha_ = fluxFunction_->gradAlpha();
}


void Foam::liquidPhaseModel::decode()
{
    this->U_ = this->alphaRhoU_/this->alphaRho_;
    this->U_.correctBoundaryConditions();

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
    thermoPtr_->correct();
}

void Foam::liquidPhaseModel::store()
{
    (*this).storeOldTime();
    alphaRho_.storeOldTime();
    alphaRhoU_.storeOldTime();
    alphaRhoE_.storeOldTime();
}

// ************************************************************************* //
