/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "YuStandishPackingLimit.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace packingLimitModels
{
    defineTypeNameAndDebug(YuStandishPackingLimit, 0);

    addToRunTimeSelectionTable
    (
        packingLimitModel,
        YuStandishPackingLimit,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::packingLimitModels::
YuStandishPackingLimit::YuStandishPackingLimit
(
    const dictionary& dict,
    const polydisperseKineticTheoryModel& kt
)
:
    packingLimitModel(dict, kt),
    residualAlpha_(dict_.lookupType<scalar>("residualAlpha"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::packingLimitModels::
YuStandishPackingLimit::~YuStandishPackingLimit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::kineticTheoryModels::packingLimitModels::
YuStandishPackingLimit::alphaMax
(
    const label celli,
    const scalarList& ds
) const
{
    scalar maxAlpha = 1.0;
    scalar alphap = kt_.alphap()[celli];

    forAll(ds, phasei)
    {
        const phaseModel& phase1 = kt_.fluid().phases()[phasei];
        scalar alpha1 = phase1[celli];
        scalar alphaMax1 = phase1.alphaMax();
        scalar d1 = ds[phasei];

        scalar cxi = alpha1/max(alphap, residualAlpha_);

        scalar sum1 = 0.0;
        scalar sum2 = 0.0;

        forAll(ds, phasej)
        {
            if (phasej == phasei)
            {
                continue;
            }
            scalar d2 = ds[phasej];

            scalar rij = d1;
            scalar Xij = 1.0;
            scalar pij = alphaMax1;

            if (phasei < phasej)
            {
                rij = d1/d2;
                Xij = (1.0 - sqr(rij))/(2.0 - alphaMax1);
            }
            else
            {
                rij = d2/d1;
                Xij = 1.0 - (1.0 - sqr(rij))/(2.0 - alphaMax1);
            }

            if (rij <= 0.741)
            {
                pij +=
                    alphaMax1
                   *(1.0 - alphaMax1)
                   *(1.0 - 2.35*rij + 1.35*sqr(rij));
            }

            scalar cmpt = (1.0 - alphaMax1/pij)*cxi/Xij;
            if (phasej < phasei)
            {
                sum1 += cmpt;
            }
            else if (phasej > phasei)
            {
                sum2 += cmpt;
            }
        }
        maxAlpha = min(maxAlpha, alphaMax1/(1.0 - sum1 - sum2));
    }
    return maxAlpha;
}


bool Foam::kineticTheoryModels::packingLimitModels::
YuStandishPackingLimit::read()
{
    return true;
}


// ************************************************************************* //
