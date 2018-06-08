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

#include "EulerPhase.H"
#include "twoPhaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFluxIntegrators
{

    defineTypeNameAndDebug(EulerPhase, 0);
    addToRunTimeSelectionTable(phaseFluxIntegrator, EulerPhase, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxIntegrators::EulerPhase::EulerPhase
(
    phaseModel& phase1,
    phaseModel& phase2
)
:
    phaseFluxIntegrator(phase1, phase2)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxIntegrators::EulerPhase::~EulerPhase()
{}


// * * * * * * * * * * * * * Public Member Fucntions * * * * * * * * * * * * //

void Foam::phaseFluxIntegrators::EulerPhase::integrateFluxes
(
    const dimensionedVector& g,
    volVectorField& Ui,
    volScalarField& pi
)
{
    dimensionedScalar deltaT = Ui.mesh().time().deltaT();

    volScalarField& alpha1 = phase1_;
    volScalarField& alpha2 = phase2_;

    phase1_.encode();
    phase2_.encode();

    phase1_.updateFluxes();
    phase2_.alphaf() = 1.0 - phase1_.alphaf();
    phase2_.updateFluxes();

    phase1_.advect
    (
        deltaT,
        g,
        Ui,
        pi,
        false
    );
    phase1_.decode();

    phase2_.advect
    (
        deltaT,
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
}