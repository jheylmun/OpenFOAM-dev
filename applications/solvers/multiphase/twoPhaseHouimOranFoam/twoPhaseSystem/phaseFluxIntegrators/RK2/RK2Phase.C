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

#include "RK2Phase.H"
#include "twoPhaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFluxIntegrators
{

    defineTypeNameAndDebug(RK2Phase, 0);
    addToRunTimeSelectionTable(phaseFluxIntegrator, RK2Phase, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxIntegrators::RK2Phase::RK2Phase
(
    phaseModel& phase1,
    phaseModel& phase2
)
:
    phaseFluxIntegrator(phase1, phase2)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxIntegrators::RK2Phase::~RK2Phase()
{}


// * * * * * * * * * * * * * Public Member Fucntions * * * * * * * * * * * * //

void Foam::phaseFluxIntegrators::RK2Phase::integrateFluxes
(
    const dimensionedVector& g,
    volVectorField& Ui,
    volScalarField& pi
)
{
    const dimensionedScalar& deltaT = Ui.mesh().time().deltaT();

    const volScalarField& alpha1 = phase1_;
    volScalarField& alpha2 = phase2_;

    //- Predictor step
    phase1_.updateFluxes();
    phase2_.alphaf() = 1.0 - phase1_.alphaf();
    phase2_.updateFluxes();

    phase1_.advect
    (
        deltaT*0.5,
        g,
        Ui,
        pi,
        true
    );
    phase1_.decode();

    phase2_.advect
    (
        deltaT*0.5,
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

    //- Corrector
    phase1_.updateFluxes();
    phase2_.alphaf() = 1.0 - phase1_.alphaf();
    phase2_.updateFluxes();

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