/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Jeff Heylmun
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

#include "HLLCPhaseFlux.H"
#include "surfaceInterpolate.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFluxFunctions
{
    defineTypeNameAndDebug(HLLCPhase, 0);
    addToRunTimeSelectionTable(phaseFluxFunction, HLLCPhase, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxFunctions::HLLCPhase::HLLCPhase
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    phaseFluxFunction(mesh, phaseName),
    residualU_("small", dimVelocity, epsilon_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxFunctions::HLLCPhase::~HLLCPhase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseFluxFunctions::HLLCPhase::updateFluxes
(
    surfaceScalarField& massFlux,
    surfaceVectorField& momentumFlux,
    surfaceScalarField& energyFlux,
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& E,
    const volScalarField& p,
    const volScalarField& a,
    const volVectorField& Ui,
    const volScalarField& pi
)
{
    surfaceVectorField normal(mesh_.Sf()/mesh_.magSf());

    surfaceScalarField alphaOwn
    (
        fvc::interpolate(alpha, own_, interpScheme(alpha.name()))
    );
    surfaceScalarField alphaNei
    (
        fvc::interpolate(alpha, nei_, interpScheme(alpha.name()))
    );

    surfaceScalarField rhoOwn
    (
        fvc::interpolate(rho, own_, interpScheme(rho.name()))
    );
    surfaceScalarField rhoNei
    (
        fvc::interpolate(rho, nei_, interpScheme(rho.name()))
    );

    surfaceVectorField UOwn(fvc::interpolate(U, own_, interpScheme(U.name())));
    surfaceVectorField UNei(fvc::interpolate(U, nei_, interpScheme(U.name())));

    surfaceScalarField EOwn(fvc::interpolate(E, own_, interpScheme(E.name())));
    surfaceScalarField ENei(fvc::interpolate(E, nei_, interpScheme(E.name())));

    surfaceScalarField pOwn(fvc::interpolate(p, own_, interpScheme(p.name())));
    surfaceScalarField pNei(fvc::interpolate(p, nei_, interpScheme(p.name())));

    surfaceScalarField UvOwn(UOwn & normal);
    surfaceScalarField UvNei(UNei & normal);

    surfaceScalarField aOwn
    (
        fvc::interpolate(a, own_, interpScheme(a.name()))
    );
    surfaceScalarField aNei
    (
        fvc::interpolate(a, nei_, interpScheme(a.name()))
    );

    // Averages
    surfaceScalarField aTilde
    (
        "aTilde",
        (sqrt(rhoOwn)*aOwn + sqrt(rhoNei)*aNei)
       /(sqrt(rhoOwn) + sqrt(rhoNei))
    );
    surfaceVectorField uTilde
    (
        "uTilde",
        (sqrt(rhoOwn)*UOwn + sqrt(rhoNei)*UNei)
       /(sqrt(rhoOwn) + sqrt(rhoNei))
    );

    // Compute wave speeds
    surfaceScalarField SOwn
    (
        "SOwn", min(UvOwn - aOwn, (uTilde & normal) - aTilde)
    );
    surfaceScalarField SNei
    (
        "SNei", max(UvNei + aNei, (uTilde & normal) + aTilde)
    );

    //- Star quantities
    surfaceScalarField SStar
    (
        "SStar",
        (
            pNei - pOwn
          + rhoOwn*UvOwn*(SOwn - UvOwn)
          - rhoNei*UvNei*(SNei - UvNei)
        )/(rhoOwn*(SOwn - UvOwn) - rhoNei*(SNei - UvNei))
    );
    surfaceScalarField pStar
    (
        "pOwnNei",
        pOwn + rhoOwn*(SOwn - UvOwn)*(SStar - UvOwn)
    );
    surfaceScalarField rhoOwnStar
    (
        "rhoOwnStar",
        rhoOwn*(SOwn - UvOwn)/(SOwn - SStar)
    );
    surfaceScalarField rhoNeiStar
    (
        "rhoNeiStar",
        rhoNei*(SNei - UvNei)/(SNei - SStar)
    );
    surfaceScalarField EOwnStar
    (
        "EOwnStar",
        EOwn + (pStar*SStar - pOwn*UvOwn)/(rhoOwn*(SOwn - UvOwn))
    );
    surfaceScalarField ENeiStar
    (
        "ENeiStar",
        ENei + (pStar*SStar - pNei*UvNei)/(rhoNei*(SNei - UvNei))
    );

    // Reimann primitive quantities
    alphaf_ =
        pos0(SOwn)*alphaOwn
      + neg(SOwn)*pos0(SStar)*alphaOwn
      + neg(SStar)*pos0(SNei)*alphaNei
      + neg(SNei)*alphaNei;

    surfaceScalarField rhoR
    (
        "rhoR",
        pos0(SOwn)*rhoOwn
      + neg(SOwn)*pos0(SStar)*rhoOwnStar
      + neg(SStar)*pos0(SNei)*rhoNeiStar
      + neg(SNei)*rhoNei
    );
    Uf_ =
        pos0(SOwn)*UOwn
      + (neg(SOwn)*pos0(SStar) + neg(SStar)*pos0(SNei))*SStar*normal
      + neg(SNei)*UNei;
    phi_ = Uf_ & mesh_.Sf();

    surfaceScalarField ER
    (
        "ER",
        pos0(SOwn)*EOwn
      + neg(SOwn)*pos0(SStar)*EOwnStar
      + neg(SStar)*pos0(SNei)*ENeiStar
      + neg(SNei)*ENei
    );
    pf_ =
        pos0(SOwn)*pOwn
      + (neg(SOwn)*pos0(SStar) + neg(SStar)*pos0(SNei))*pStar
      + neg(SNei)*pNei;

    // Set total fluxes
    massFlux = alphaf_*rhoR*phi();
    momentumFlux = massFlux*Uf_ + alphaf_*pf_*mesh_.Sf();
    energyFlux = alphaf_*phi()*(rhoR*ER + pf_);
}


void Foam::phaseFluxFunctions::HLLCPhase::updateFluxes
(
    surfaceScalarField& massFlux,
    surfaceVectorField& momentumFlux,
    surfaceScalarField& energyFlux,
    const surfaceScalarField& alphaf,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& E,
    const volScalarField& p,
    const volScalarField& a,
    const volVectorField& Ui,
    const volScalarField& pi
)
{
    surfaceVectorField normal(mesh_.Sf()/mesh_.magSf());
    alphaf_ = alphaf;

    surfaceScalarField rhoOwn
    (
        fvc::interpolate(rho, own_, interpScheme(rho.name()))
    );
    surfaceScalarField rhoNei
    (
        fvc::interpolate(rho, nei_, interpScheme(rho.name()))
    );

    surfaceVectorField UOwn(fvc::interpolate(U, own_, interpScheme(U.name())));
    surfaceVectorField UNei(fvc::interpolate(U, nei_, interpScheme(U.name())));

    surfaceScalarField EOwn(fvc::interpolate(E, own_, interpScheme(E.name())));
    surfaceScalarField ENei(fvc::interpolate(E, nei_, interpScheme(E.name())));

    surfaceScalarField pOwn(fvc::interpolate(p, own_, interpScheme(p.name())));
    surfaceScalarField pNei(fvc::interpolate(p, nei_, interpScheme(p.name())));

    surfaceScalarField UvOwn(UOwn & normal);
    surfaceScalarField UvNei(UNei & normal);

    surfaceScalarField aOwn(fvc::interpolate(a, own_, interpScheme(a.name())));
    surfaceScalarField aNei(fvc::interpolate(a, nei_, interpScheme(a.name())));

    // Averages
    surfaceScalarField aTilde
    (
        "aTilde",
        (sqrt(rhoOwn)*aOwn + sqrt(rhoNei)*aNei)
       /(sqrt(rhoOwn) + sqrt(rhoNei))
    );
    surfaceScalarField uTilde
    (
        "uTilde",
        (
            (sqrt(rhoOwn)*UOwn + sqrt(rhoNei)*UNei)
           /(sqrt(rhoOwn) + sqrt(rhoNei))
        ) & normal
    );

    // Compute wave speeds
    surfaceScalarField SOwn
    (
        "SOwn", min(UvOwn - aOwn, uTilde - aTilde)
    );
    surfaceScalarField SNei
    (
        "SNei", max(UvNei + aNei, uTilde + aTilde)
    );

    //- Star quantities
    surfaceScalarField SStar
    (
        "SStar",
        (
            pNei - pOwn
          + rhoOwn*UvOwn*(SOwn - UvOwn)
          - rhoNei*UvNei*(SNei - UvNei)
        )/(rhoOwn*(SOwn - UvOwn) - rhoNei*(SNei - UvNei))
    );
    surfaceScalarField pStar
    (
        "pOwnNei",
        (
            pOwn + rhoOwn*(SOwn - UvOwn)*(SStar - UvOwn)
          + pNei + rhoNei*(SNei - UvNei)*(SStar - UvNei)
        )*0.5
    );
    surfaceScalarField rhoOwnStar
    (
        "rhoOwnStar",
        rhoOwn*(SOwn - UvOwn)/(SOwn - SStar)
    );
    surfaceScalarField rhoNeiStar
    (
        "rhoNeiStar",
        rhoNei*(SNei - UvNei)/(SNei - SStar)
    );
    surfaceScalarField EOwnStar
    (
        "EOwnStar",
        EOwn + (pStar*SStar - pOwn*UvOwn)/(rhoOwn*(SOwn - UvOwn))
    );
    surfaceScalarField ENeiStar
    (
        "ENeiStar",
        ENei + (pStar*SStar - pNei*UvNei)/(rhoNei*(SNei - UvNei))
    );

    // Reimann primitive quantities
    surfaceScalarField rhoR
    (
        "rhoR",
        pos0(SOwn)*rhoOwn
      + neg(SOwn)*pos0(SStar)*rhoOwnStar
      + neg(SStar)*pos0(SNei)*rhoNeiStar
      + neg(SNei)*rhoNei
    );
    Uf_ =
        pos0(SOwn)*UOwn
      + (neg(SOwn)*pos0(SStar) + neg(SStar)*pos0(SNei))*SStar*normal
      + neg(SNei)*UNei;
    phi_ = Uf_ & mesh_.Sf();

    surfaceScalarField ER
    (
        "ER",
        pos0(SOwn)*EOwn
      + neg(SOwn)*pos0(SStar)*EOwnStar
      + neg(SStar)*pos0(SNei)*ENeiStar
      + neg(SNei)*ENei
    );

    pf_ =
        pos0(SOwn)*pOwn
      + (neg(SOwn)*pos0(SStar) + neg(SStar)*pos0(SNei))*pStar
      + neg(SNei)*pNei;

    // Set total fluxes

    massFlux = alphaf_*rhoR*phi();
    momentumFlux = massFlux*Uf_ + alphaf_*pf_*mesh_.Sf();
    energyFlux = alphaf_*phi()*(rhoR*ER + pf_);
}