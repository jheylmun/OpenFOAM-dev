/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "segregated.H"
#include "phasePair.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(segregated, 0);
    addToRunTimeSelectionTable(dragModel, segregated, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::segregated::segregated
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    m_("m", dimless, dict),
    n_("n", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::segregated::~segregated()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::segregated::CdRe() const
{
    FatalErrorInFunction
        << "Not implemented."
        << "Drag coefficient not defined for the segregated model."
        << exit(FatalError);

    return pair_.phase1();
}


Foam::tmp<Foam::volScalarField> Foam::dragModels::segregated::K() const
{
    const fvMesh& mesh(pair_.phase1().mesh());

    const volScalarField& alpha1(pair_.phase1());
    const volScalarField& alpha2(pair_.phase2());

    const volScalarField& rho1(pair_.phase1().rho());
    const volScalarField& rho2(pair_.phase2().rho());

    tmp<volScalarField> tnu1(pair_.phase1().nu());
    tmp<volScalarField> tnu2(pair_.phase2().nu());

    const volScalarField& nu1(tnu1());
    const volScalarField& nu2(tnu2());

    volScalarField L
    (
        IOobject
        (
            "L",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("L", dimLength, 0),
        zeroGradientFvPatchField<scalar>::typeName
    );
    L.primitiveFieldRef() = cbrt(mesh.V());
    L.correctBoundaryConditions();

    volScalarField I
    (
        alpha1
       /max
        (
            alpha1 + alpha2,
            pair_.phase1().residualAlpha() + pair_.phase2().residualAlpha()
        )
    );
    volScalarField magGradI
    (
        max
        (
            mag(fvc::grad(I)),
            (pair_.phase1().residualAlpha() + pair_.phase2().residualAlpha())/L
        )
    );

    volScalarField muI
    (
        rho1*nu1*rho2*nu2
       /(rho1*nu1 + rho2*nu2)
    );
    volScalarField muAlphaI
    (
        alpha1*rho1*nu1*alpha2*rho2*nu2
       /(
           max(alpha1, pair_.phase1().residualAlpha())*rho1*nu1
         + max(alpha2, pair_.phase2().residualAlpha())*rho2*nu2
        )
    );

    volScalarField ReI
    (
        pair_.rho()
       *pair_.magUr()
       /(magGradI*muI)
    );

    volScalarField lambda(m_*ReI + n_*muAlphaI/muI);

    return lambda*sqr(magGradI)*muI;
}


Foam::tmp<Foam::surfaceScalarField> Foam::dragModels::segregated::Kf() const
{
    return fvc::interpolate(K());
}


Foam::scalar Foam::dragModels::segregated::cellCdRe
(
    const label celli
) const
{
    FatalErrorInFunction
        << "Not implemented."
        << "Drag coefficient not defined for the segregated model."
        << exit(FatalError);

    return 0.0;
}


Foam::scalar Foam::dragModels::segregated::cellK
(
    const label celli
) const
{
    const fvMesh& mesh(pair_.phase1().mesh());

    scalar alpha1(pair_.phase1()[celli]);
    scalar alpha2(pair_.phase2()[celli]);

    scalar rho1(pair_.phase1().thermo().cellrho(celli));
    scalar rho2(pair_.phase2().thermo().cellrho(celli));

    scalar nu1(pair_.phase1().thermo().cellnu(celli));
    scalar nu2(pair_.phase2().thermo().cellnu(celli));

    vector gradI(Zero);
    const labelList& cell = mesh.cells()[celli];
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const label nInternalFaces = mesh.nInternalFaces();
    forAll(cell, facei)
    {
        if (cell[facei] > nInternalFaces) continue;
        label own = owner[cell[facei]];
        label nei = neighbour[cell[facei]];
        scalar alpha1Own(pair_.phase1()[own]);
        scalar alpha1Nei(pair_.phase1()[nei]);
        scalar alpha2Own(pair_.phase2()[own]);
        scalar alpha2Nei(pair_.phase2()[nei]);

        scalar IOwn
        (
            alpha1Own
           /max
            (
                alpha1Own + alpha2Own,
                pair_.phase1().residualAlpha().value()
              + pair_.phase2().residualAlpha().value()
            )
        );
        scalar INei
        (
            alpha1Nei
           /max
            (
                alpha1Nei + alpha2Nei,
                pair_.phase1().residualAlpha().value()
              + pair_.phase2().residualAlpha().value()
            )
        );
        gradI += (IOwn - INei)*mesh.Sf()[cell[facei]];
    }
    scalar L = cbrt(mesh.V()[celli]);

    scalar magGradI
    (
        max
        (
            mag(gradI),
            (
                pair_.phase1().residualAlpha().value()
              + pair_.phase2().residualAlpha().value()
            )/L
        )
    );

    scalar muI(rho1*nu1*rho2*nu2/(rho1*nu1 + rho2*nu2));

    scalar muAlphaI
    (
        alpha1*rho1*nu1*alpha2*rho2*nu2
       /(
            max
            (
                alpha1,
                pair_.phase1().residualAlpha().value()
            )*rho1*nu1
          + max
            (
                alpha2,
                pair_.phase2().residualAlpha().value()
            )*rho2*nu2
        )
    );

    scalar ReI
    (
        pair_.rho(celli)
       *pair_.magUr(celli)
       /(magGradI*muI)
    );

    scalar lambda(m_.value()*ReI + n_.value()*muAlphaI/muI);

    return lambda*sqr(magGradI)*muI;
}


// ************************************************************************* //
