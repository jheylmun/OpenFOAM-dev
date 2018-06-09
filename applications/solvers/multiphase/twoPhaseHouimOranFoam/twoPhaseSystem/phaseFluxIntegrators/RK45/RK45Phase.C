/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "RK45Phase.H"
#include "twoPhaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFluxIntegrators
{

    defineTypeNameAndDebug(RK45Phase, 0);
    addToRunTimeSelectionTable(phaseFluxIntegrator, RK45Phase, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxIntegrators::RK45Phase::RK45Phase
(
    phaseModel& phase1,
    phaseModel& phase2
)
:
    phaseFluxIntegrator(phase1, phase2)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxIntegrators::RK45Phase::~RK45Phase()
{}


// * * * * * * * * * * * * * Public Member Fucntions * * * * * * * * * * * * //

void Foam::phaseFluxIntegrators::RK45Phase::integrateFluxes
(
    const dimensionedVector& g,
    volVectorField& Ui,
    volScalarField& pi
)
{
    const dimensionedScalar& deltaT = Ui.mesh().time().deltaT();
    const volScalarField& alpha1 = phase1_;
    volScalarField& alpha2 = phase2_;

    //- 1st predictor step
    phase1_.updateFluxes();
    surfaceScalarField massFlux1(phase1_.massFlux());
    surfaceVectorField momentumFlux1(phase1_.momentumFlux());
    surfaceScalarField energyFlux1(phase1_.energyFlux());
    tmp<surfaceScalarField> PTEFlux1;
    if (phase1_.granular())
    {
        PTEFlux1 = tmp<surfaceScalarField>
        (
            new surfaceScalarField(phase1_.PTEFlux())
        );
    }

    phase2_.alphaf() = 1.0 - phase1_.alphaf();
    phase2_.updateFluxes();
    surfaceScalarField massFlux2(phase2_.massFlux());
    surfaceVectorField momentumFlux2(phase2_.momentumFlux());
    surfaceScalarField energyFlux2(phase2_.energyFlux());
    tmp<surfaceScalarField> PTEFlux2;
    if (phase2_.granular())
    {
        PTEFlux2 = tmp<surfaceScalarField>
        (
            new surfaceScalarField(phase2_.PTEFlux())
        );
    }

    volVectorField KUi(Ui);
    volScalarField Kpi(pi);
    volVectorField KgradAlpha(gradAlpha_);

    phase1_.advect
    (
        deltaT*0.5,
        g,
        Ui,
        pi,
        false
    );
    phase1_.decode();

    phase2_.advect
    (
        deltaT*0.5,
        g,
        Ui,
        pi,
        false
    );
    alpha2 = 1.0 - alpha1;
    alpha2.correctBoundaryConditions();
    phase2_.decode();


    Ui = phase1_.fluid().mixtureU();
    pi = phase1_.fluid().mixturep();

    //- 2nd predictor step
    phase1_.updateFluxes();
    massFlux1 += 2.0*phase1_.massFlux();
    momentumFlux1 += 2.0*phase1_.momentumFlux();
    energyFlux1 += 2.0*phase1_.energyFlux();
    if (phase1_.granular())
    {
        PTEFlux1.ref() += 2.0*phase1_.PTEFlux();
    }

    phase2_.alphaf() = 1.0 - phase1_.alphaf();
    phase2_.updateFluxes();
    massFlux2 += 2.0*phase2_.massFlux();
    momentumFlux2 += 2.0*phase2_.momentumFlux();
    energyFlux2 += 2.0*phase2_.energyFlux();
    if (phase2_.granular())
    {
        PTEFlux2.ref() += 2.0*phase2_.PTEFlux();
    }

    KUi += 2.0*Ui;
    Kpi += 2.0*pi;
    KgradAlpha += 2.0*gradAlpha_;

    phase1_.advect
    (
        deltaT*0.5,
        g,
        Ui,
        pi,
        false
    );
    phase1_.decode();

    phase2_.advect
    (
        deltaT*0.5,
        g,
        Ui,
        pi,
        false
    );
    alpha2 = 1.0 - alpha1;
    alpha2.correctBoundaryConditions();
    phase2_.decode();

    Ui = phase1_.fluid().mixtureU();
    pi = phase1_.fluid().mixturep();

    //- Third predictor step
    phase1_.updateFluxes();
    massFlux1 += 2.0*phase1_.massFlux();
    momentumFlux1 += 2.0*phase1_.momentumFlux();
    energyFlux1 += 2.0*phase1_.energyFlux();
    if (phase1_.granular())
    {
        PTEFlux1.ref() += 2.0*phase1_.PTEFlux();
    }

    phase2_.alphaf() = 1.0 - phase1_.alphaf();
    phase2_.updateFluxes();
    massFlux2 += 2.0*phase2_.massFlux();
    momentumFlux2 += 2.0*phase2_.momentumFlux();
    energyFlux2 += 2.0*phase2_.energyFlux();
    if (phase2_.granular())
    {
        PTEFlux2.ref() += 2.0*phase2_.PTEFlux();
    }

    KUi += 2.0*Ui;
    Kpi += 2.0*pi;
    KgradAlpha += 2.0*gradAlpha_;

    phase1_.advect
    (
        deltaT,
        g,
        Ui,
        pi,
        true
    );
    phase1_.decode();

    phase2_.advect
    (
        deltaT,
        g,
        Ui,
        pi,
        true
    );
    alpha2 = 1.0 - alpha1;
    alpha2.correctBoundaryConditions();
    phase2_.decode();

    Ui = phase1_.fluid().mixtureU();
    pi = phase1_.fluid().mixturep();

    //- Final correction step
    phase1_.updateFluxes();
    phase1_.massFlux() = (massFlux1 + phase1_.massFlux())/6.0;
    phase1_.momentumFlux() = (momentumFlux1 + phase1_.momentumFlux())/6.0;
    phase1_.energyFlux() = (energyFlux1 + phase1_.energyFlux())/6.0;
    if (phase1_.granular())
    {
        phase1_.PTEFlux() = (PTEFlux1 + phase1_.PTEFlux())/6.0;
    }

    phase2_.alphaf() = 1.0 - phase1_.alphaf();
    phase2_.updateFluxes();
    phase2_.massFlux() = (massFlux2 + phase2_.massFlux())/6.0;
    phase2_.momentumFlux() = (momentumFlux2 + phase2_.momentumFlux())/6.0;
    phase2_.energyFlux() = (energyFlux2 + phase2_.energyFlux())/6.0;
    if (phase2_.granular())
    {
        phase2_.PTEFlux() = (PTEFlux2 + phase2_.PTEFlux())/6.0;
    }

    Ui = (KUi + Ui)/6.0;
    pi = (Kpi + pi)/6.0;
    phase1_.gradAlpha() = (KgradAlpha + gradAlpha_)/6.0;
    phase2_.gradAlpha() = -gradAlpha_;

    phase1_.advect
    (
        deltaT,
        g,
        Ui,
        pi,
        true
    );
    phase1_.decode();

    phase2_.advect
    (
        deltaT,
        g,
        Ui,
        pi,
        true
    );
    alpha2 = 1.0 - alpha1;
    alpha2.correctBoundaryConditions();
    phase2_.decode();

    Ui = phase1_.fluid().mixtureU();
    pi = phase1_.fluid().mixturep();
}