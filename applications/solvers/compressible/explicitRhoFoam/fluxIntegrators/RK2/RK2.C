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

#include "RK2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxIntegrators
{

    defineTypeNameAndDebug(RK2, 0);
    addToRunTimeSelectionTable(fluxIntegrator, RK2, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxIntegrators::RK2::RK2
(
    compressibleSystem& fluid,
    const word& phaseName
)
:
    fluxIntegrator(fluid, phaseName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxIntegrators::RK2::~RK2()
{}


// * * * * * * * * * * * * * Public Member Fucntions * * * * * * * * * * * * //

void Foam::fluxIntegrators::RK2::integrateFluxes
(
    volScalarField& rho,
    volVectorField& rhoU,
    volScalarField& rhoE
)
{
    fluid_.encode();
    const dimensionedScalar& deltaT = rho.mesh().time().deltaT();

    //- Predictor step
    fluid_.updateFluxes();

    rho -=  0.5*deltaT*fvc::div(massFlux_);
    rhoU -= 0.5*deltaT*fvc::div(momentumFlux_);
    rhoE -= 0.5*deltaT*fvc::div(energyFlux_);
    rho.correctBoundaryConditions();
    rhoU.correctBoundaryConditions();
    rhoE.correctBoundaryConditions();
    fluid_.decode();

    //- Corrector
    fluid_.updateFluxes();
}