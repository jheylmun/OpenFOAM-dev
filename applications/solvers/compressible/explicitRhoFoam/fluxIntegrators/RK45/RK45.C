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

#include "RK45.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxIntegrators
{

    defineTypeNameAndDebug(RK45, 0);
    addToRunTimeSelectionTable(fluxIntegrator, RK45, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxIntegrators::RK45::RK45
(
    compressibleSystem& fluid,
    const word& phaseName
)
:
    fluxIntegrator(fluid, phaseName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxIntegrators::RK45::~RK45()
{}


// * * * * * * * * * * * * * Public Member Fucntions * * * * * * * * * * * * //

void Foam::fluxIntegrators::RK45::integrateFluxes
(
    volScalarField& rho,
    volVectorField& rhoU,
    volScalarField& rhoE
)
{
    fluid_.encode();
    const dimensionedScalar& deltaT = rho.mesh().time().deltaT();

    volScalarField rhoOld(rho);
    volVectorField rhoUOld(rhoU);
    volScalarField rhoEOld(rhoE);

    //- 1st predictor step
    fluid_.updateFluxes();
    surfaceScalarField K1rhoFlux = massFlux_;
    surfaceVectorField K1rhoUFlux = momentumFlux_;
    surfaceScalarField K1rhoEFlux = energyFlux_;

    rho -=  0.5*deltaT*fvc::div(K1rhoFlux);
    rhoU -= 0.5*deltaT*fvc::div(K1rhoUFlux);
    rhoE -= 0.5*deltaT*fvc::div(K1rhoEFlux);
    rho.correctBoundaryConditions();
    rhoU.correctBoundaryConditions();
    rhoE.correctBoundaryConditions();
    fluid_.decode();

    //- 2nd predictor step
    fluid_.updateFluxes();
    surfaceScalarField K2rhoFlux = massFlux_;
    surfaceVectorField K2rhoUFlux = momentumFlux_;
    surfaceScalarField K2rhoEFlux = energyFlux_;

    rho -= 0.5*deltaT*fvc::div(K2rhoFlux);
    rhoU -= 0.5*deltaT*fvc::div(K2rhoUFlux);
    rhoE -= 0.5*deltaT*fvc::div(K2rhoEFlux);
    rho.correctBoundaryConditions();
    rhoU.correctBoundaryConditions();
    rhoE.correctBoundaryConditions();
    fluid_.decode();

    //- Third predictor step
    fluid_.updateFluxes();
    surfaceScalarField K3rhoFlux = massFlux_;
    surfaceVectorField K3rhoUFlux = momentumFlux_;
    surfaceScalarField K3rhoEFlux = energyFlux_;

    rho = rhoOld - deltaT*fvc::div(K3rhoFlux);
    rhoU = rhoUOld - deltaT*fvc::div(K3rhoUFlux);
    rhoE = rhoEOld - deltaT*fvc::div(K3rhoEFlux);
    rho.correctBoundaryConditions();
    rhoU.correctBoundaryConditions();
    rhoE.correctBoundaryConditions();
    fluid_.decode();

    //- Final correction step
    fluid_.updateFluxes();
    massFlux_ =
    (
        K1rhoFlux
      + K2rhoFlux*2.0
      + K3rhoFlux*2.0
      + massFlux_
    )/6.0;

    momentumFlux_ =
    (
        K1rhoUFlux
      + K2rhoUFlux*2.0
      + K3rhoUFlux*2.0
      + momentumFlux_
    )/6.0;

    energyFlux_ =
    (
        K1rhoEFlux
      + K2rhoEFlux*2.0
      + K3rhoEFlux*2.0
      + energyFlux_
    )/6.0;
}
