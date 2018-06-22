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
    if (!otherPhase().granular())
    {
        alphas_[stepi] = *this;
        deltaAlpha_[stepi] = -(Ui & gradAlpha_);
    }
    alphaRhos_[stepi] = alphaRho_;
    alphaRhoUs_[stepi] = alphaRhoU_;
    alphaRhoEs_[stepi] = alphaRhoE_;

    deltaAlphaRho_[stepi] = -fvc::div(massFlux_);
    deltaAlphaRhoU_[stepi] =
       - fvc::div(momentumFlux_)
       + pi*gradAlpha_
       + F
       + alphaRho_*g;
    deltaAlphaRhoE_[stepi] = -fvc::div(energyFlux_) + ((*this)*rho_*g & U_);

    if (otherPhase().granular())
    {
        deltaAlphaRhoE_[stepi] +=
            (F & otherPhase().U())
          - pi*fvc::div(otherPhase().alphaPhi());
    }
    tmp<volScalarField> alpha;
    tmp<volScalarField> deltaAlpha;
    if (!otherPhase().granular())
    {
        alpha.ref() = (*this)*coeffs[stepi];
        deltaAlpha.ref() = deltaAlpha_[stepi];
    }
    volScalarField alphaRho(alphaRho_*coeffs[stepi]);
    volVectorField alphaRhoU(alphaRhoU_*coeffs[stepi]);
    volScalarField alphaRhoE(alphaRhoE_*coeffs[stepi]);

    volScalarField deltaAlphaRho(deltaAlphaRho_[stepi]*Fcoeffs[stepi]);
    volVectorField deltaAlphaRhoU(deltaAlphaRhoU_[stepi]*Fcoeffs[stepi]);
    volScalarField deltaAlphaRhoE(deltaAlphaRhoE_[stepi]*Fcoeffs[stepi]);

    for (label i = 0; i < stepi; i++)
    {
        if (!otherPhase().granular())
        {
            alpha.ref() += alphas_[i]*coeffs[i];
            deltaAlpha.ref() += deltaAlpha_[i]*Fcoeffs[i];
        }
        alphaRho += alphaRhos_[i]*coeffs[i];
        alphaRhoU += alphaRhoUs_[i]*coeffs[i];
        alphaRhoE += alphaRhoEs_[i]*coeffs[i];

        deltaAlphaRho += deltaAlphaRho_[i]*Fcoeffs[i];
        deltaAlphaRhoU += deltaAlphaRhoU_[i]*Fcoeffs[i];
        deltaAlphaRhoE += deltaAlphaRhoE_[i]*Fcoeffs[i];
    }
    if (!otherPhase().granular())
    {
        refCast<volScalarField>(*this) += deltaT*deltaAlpha;
        this->correctBoundaryConditions();
    }

    alphaRho_ = alphaRho + deltaT*deltaAlphaRho;
    alphaRho_.correctBoundaryConditions();

    alphaRhoU_ = alphaRhoU + deltaT*deltaAlphaRhoU;
    alphaRhoU_.correctBoundaryConditions();

    alphaRhoE_ = alphaRhoE + deltaT*deltaAlphaRhoE;
    alphaRhoE_.correctBoundaryConditions();
}


void Foam::fluidPhaseModel::updateFluxes()
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
        E_,
        p_,
        c(),
        fluid_.U(),
        fluid_.p()
    );
    gradAlpha_ = fluxFunction_->gradAlpha();
}

void Foam::fluidPhaseModel::updateFluxes(const surfaceScalarField& alphaf)
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
        E_,
        p_,
        c(),
        fluid_.U(),
        fluid_.p()
    );
    gradAlpha_ = fluxFunction_->gradAlpha();
}


void Foam::fluidPhaseModel::decode()
{
    if (otherPhase().slavePressure())
    {
        const dictionary& eosDict =
            fluid_.mesh().lookupObject<dictionary>
            (
                IOobject::groupName
                (
                    "thermophysicalProperties",
                    otherPhase().name()
                )
            ).subDict("mixture").subDict("equationOfState");
        dimensionedScalar pInf
        (
            "p0",
            dimPressure,
            eosDict
        );
        dimensionedScalar otherGamma
        (
            "gamma",
            dimless,
            eosDict
        );

        U_ = alphaRhoU_/alphaRho_;
        U_.correctBoundaryConditions();

        E_ = alphaRhoE_/alphaRho_;
        he_ = E_ - 0.5*magSqr(U_);
        he_.correctBoundaryConditions();

        phaseModel& liquid = fluid_.mesh().lookupObjectRef<phaseModel>
        (
            IOobject::groupName("alpha", otherPhase().name())
        );
        volScalarField& alpha2 = liquid;
        const rhoThermo& otherThermo = liquid.thermo();

        volScalarField A(alphaRho_*(thermoPtr_->gamma() - 1.0)*he_);
        volScalarField B
        (
            liquid.alphaRho()*(otherGamma - 1.0)*otherThermo.he()
        );
        volScalarField delta
        (
            sqr(otherGamma*pInf - A - B) + 4.0*A*otherGamma*pInf
        );
        thermoPtr_->p() = 0.5*(A + B - otherGamma*pInf + sqrt(delta));
        thermoPtr_->correct();

        p_ = thermoPtr_->p();
        p_.correctBoundaryConditions();

        volScalarField& alpha = *this;
        alpha = A/thermoPtr_->p();
        alpha.max(residualAlpha_);
        rho_ = alphaRho_/alpha;

        alpha2 = 1.0 - alpha;
        liquid.rho() =
            liquid.alphaRho()/Foam::max(alpha2, liquid.residualAlpha());
        liquid.thermo().correct();
    }
    else
    {
        this->max(residualAlpha_);
        rho_ = alphaRho_/(*this);
        rho_.correctBoundaryConditions();

        U_ = alphaRhoU_/alphaRho_;
        U_.correctBoundaryConditions();

        E_ = alphaRhoE_/alphaRho_;
        E_.correctBoundaryConditions();

        he_ = E_ - 0.5*magSqr(U_);
        he_.correctBoundaryConditions();

        thermoPtr_->correctP();
        p_ = thermoPtr_->p();
        p_.correctBoundaryConditions();
    }
}


void Foam::fluidPhaseModel::encode()
{
    (*this).max(residualAlpha_);

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
    E_.correctBoundaryConditions();

    alphaRhoE_ = alphaRho_*E_;
    alphaRhoE_.correctBoundaryConditions();

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
