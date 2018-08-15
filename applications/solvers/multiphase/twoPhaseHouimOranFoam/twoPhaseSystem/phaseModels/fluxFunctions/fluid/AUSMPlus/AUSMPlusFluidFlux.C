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

#include "AUSMPlusFluidFlux.H"
#include "surfaceInterpolate.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidFluxFunctions
{
    defineTypeNameAndDebug(AUSMPlusFlux, 0);
    addToRunTimeSelectionTable(fluidFluxFunction, AUSMPlusFlux, dictionary);
}
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::fluidFluxFunctions::AUSMPlusFlux::M1
(
    const surfaceScalarField& M,
    const label sign
)
{
    return 0.5*(M + sign*mag(M));
}

Foam::tmp<Foam::surfaceScalarField>
Foam::fluidFluxFunctions::AUSMPlusFlux::M2
(
    const surfaceScalarField& M,
    const label sign
)
{
    return sign*0.25*sqr(M + sign);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidFluxFunctions::AUSMPlusFlux::AUSMPlusFlux
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    fluidFluxFunction(mesh, phaseName),
    fa_(dict_.lookupOrDefault("fa", 1.0)),
    D_(dict_.lookupOrDefault("D", 1.0)),
    alphaMax_(dict_.lookupOrDefault("alphaMax", 0.63)),
    alphaMinFriction_(dict_.lookupOrDefault("alphaMinFriction", 0.5)),
    cutOffMa_("small", dimless, epsilon_),
    residualRho_("small", dimDensity, epsilon_),
    residualU_("small", dimVelocity, epsilon_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidFluxFunctions::AUSMPlusFlux::~AUSMPlusFlux()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluidFluxFunctions::AUSMPlusFlux::updateFluxes
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

    surfaceScalarField HOwn(EOwn + pOwn/rhoOwn);
    surfaceScalarField HNei(ENei + pNei/rhoNei);

    surfaceScalarField aOwn(fvc::interpolate(a, own_, interpScheme(a.name())));
    surfaceScalarField aNei(fvc::interpolate(a, nei_, interpScheme(a.name())));

    surfaceScalarField UvOwn(UOwn & normal);
    surfaceScalarField UvNei(UNei & normal);

    surfaceScalarField aStar(0.5*(aOwn + aNei));

    // Compute slpit Mach numbers
    surfaceScalarField MaOwn("MaOwn", UvOwn/aStar);
    surfaceScalarField MaNei("MaNei", UvNei/aStar);
    surfaceScalarField magMaOwn(mag(MaOwn));
    surfaceScalarField magMaNei(mag(MaNei));

    surfaceScalarField M0
    (
        "M0",
        sqrt
        (
            min
            (
                1.0,
                sqr((MaOwn + MaNei)*0.5)
            )
        )
    );
    surfaceScalarField fa(M0*(2.0 - M0));
    surfaceScalarField A(3.0/16.0*(5.0*sqr(fa) - 4.0));

    surfaceScalarField BetaOwn
    (
        "BetaOwn",
        pos0(magMaOwn - 1)*0.5*(1.0 + sign(MaOwn))
      + neg(magMaOwn - 1)
       *(
            0.25*(2.0 - MaOwn)*sqr(MaOwn + 1.0)
          + A*MaOwn*sqr(sqr(MaOwn) - 1.0)
        )
    );

    surfaceScalarField BetaNei
    (
        "BetaNei",
        pos0(magMaNei - 1)*0.5*(1.0 - sign(MaNei))
      + neg(magMaNei - 1)
       *(
            0.25*(2.0 + MaNei)*sqr(MaNei - 1.0)
          - A*MaNei*sqr(sqr(MaNei) - 1.0)
        )
    );

    surfaceScalarField magUBar
    (
        (alphaOwn*rhoOwn*mag(UvOwn) + alphaNei*rhoNei*mag(UvNei))
       /(alphaOwn*rhoOwn + alphaNei*rhoNei)
    );

    surfaceScalarField Mhat
    (
        min(1.0, sqrt((sqr(UvOwn) + sqr(UvNei))*0.5)/aStar)
    );
    surfaceScalarField Xi(sqr(1.0 + Mhat));
    surfaceScalarField gOwn(-max(min(MaOwn, 0.0), -1.0));
    surfaceScalarField gNei(min(max(MaNei, 0.0), 1.0));

    surfaceScalarField mDot
    (
        "mDot",
        0.5*(alphaOwn*rhoOwn*UvOwn + alphaNei*rhoNei*UvNei)
      - 0.5*magUBar
       *(
            alphaNei*rhoNei*(1.0 - gNei)
          - alphaOwn*rhoOwn*(1.0 - gOwn)
        )
      - Xi*(alphaNei*pNei - alphaOwn*pOwn)/(2.0*aStar)
    );

    alphaf_ = pos0(mDot)*alphaOwn + neg(mDot)*alphaNei;
    pf_ =
        (BetaOwn*alphaOwn*pOwn + BetaNei*alphaNei*pNei)
       /max(alphaOwn + alphaNei, residualAlpha_);

    massFlux = mesh_.magSf()*mDot;
    Uf_ = mDot/(0.5*(rhoOwn + rhoNei))*normal;
    phi_ = Uf_ & mesh_.Sf();

    momentumFlux =
        mesh_.magSf()*0.5
       *(
            mDot*(UOwn + UNei)
          + mag(mDot)*(UOwn - UNei)
        )
      + (
            BetaOwn*alphaOwn*pOwn
          + BetaNei*alphaNei*pNei
        )*mesh_.Sf();

    energyFlux =
        mesh_.magSf()*0.5
       *(
            mDot*(HOwn + HNei)
          + mag(mDot)*(HOwn - HNei)
        );
}


void Foam::fluidFluxFunctions::AUSMPlusFlux::updateFluxes
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

    surfaceScalarField HOwn(EOwn + pOwn/rhoOwn);
    surfaceScalarField HNei(ENei + pNei/rhoNei);

    surfaceScalarField aOwn(fvc::interpolate(a, own_, interpScheme(a.name())));
    surfaceScalarField aNei(fvc::interpolate(a, nei_, interpScheme(a.name())));

    surfaceScalarField UvOwn(UOwn & normal);
    surfaceScalarField UvNei(UNei & normal);

    surfaceScalarField aStar(0.5*(aOwn + aNei));

    // Compute slpit Mach numbers
    surfaceScalarField MaOwn("MaOwn", UvOwn/aStar);
    surfaceScalarField MaNei("MaNei", UvNei/aStar);
    surfaceScalarField magMaOwn(mag(MaOwn));
    surfaceScalarField magMaNei(mag(MaNei));

    surfaceScalarField M0
    (
        "M0",
        sqrt(min(1.0, sqr((MaOwn + MaNei)*0.5)))
    );
    surfaceScalarField fa(M0*(2.0 - M0));
    surfaceScalarField A(3.0/16.0*(5.0*sqr(fa) - 4.0));

    surfaceScalarField BetaOwn
    (
        "BetaOwn",
        pos0(magMaOwn - 1)*0.5*(1.0 + sign(MaOwn))
      + neg(magMaOwn - 1)
       *(
            0.25*(2.0 - MaOwn)*sqr(MaOwn + 1.0)
          + A*MaOwn*sqr(sqr(MaOwn) - 1.0)
        )
    );

    surfaceScalarField BetaNei
    (
        "BetaNei",
        pos0(magMaNei - 1)*0.5*(1.0 - sign(MaNei))
      + neg(magMaNei - 1)
       *(
            0.25*(2.0 + MaNei)*sqr(MaNei - 1.0)
          - A*MaNei*sqr(sqr(MaNei) - 1.0)
        )
    );

    surfaceScalarField magUBar
    (
        (rhoOwn*mag(UvOwn) + rhoNei*mag(UvNei))/(rhoOwn + rhoNei)
    );

    surfaceScalarField Mhat
    (
        min(1.0, sqrt((sqr(UvOwn) + sqr(UvNei))*0.5)/aStar)
    );
    surfaceScalarField Xi(sqr(1.0 + Mhat));
    surfaceScalarField gOwn(-max(min(MaOwn, 0.0), -1.0));
    surfaceScalarField gNei(min(max(MaNei, 0.0), 1.0));

    surfaceScalarField mDot
    (
        "mDot",
        (
            0.5*(rhoOwn*UvOwn + rhoNei*UvNei)
          - 0.5*magUBar*(rhoNei*(1.0 - gNei) - rhoOwn*(1.0 - gOwn))
          - Xi*(pNei - pOwn)/(2.0*aStar)
        )*alphaf
    );

    pf_ = (BetaOwn*pOwn + BetaNei*pNei);
    massFlux = mesh_.magSf()*mDot;
    Uf_ = mDot/(0.5*(rhoOwn + rhoNei))*normal;
    phi_ = Uf_ & mesh_.Sf();

    momentumFlux =
        mesh_.magSf()*0.5
       *(
            mDot*(UOwn + UNei)
          + mag(mDot)*(UOwn - UNei)
        )
      + (
            BetaOwn*pOwn + BetaNei*pNei
        )*alphaf*mesh_.Sf();

    energyFlux =
        mesh_.magSf()*0.5
       *(
            mDot*(HOwn + HNei)
          + mag(mDot)*(HOwn - HNei)
        );
}