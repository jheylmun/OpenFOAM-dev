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

#include "bidispersePackingLimit.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace packingLimitModels
{
    defineTypeNameAndDebug(bidispersePackingLimit, 0);

    addToRunTimeSelectionTable
    (
        packingLimitModel,
        bidispersePackingLimit,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::packingLimitModels::
bidispersePackingLimit::bidispersePackingLimit
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
bidispersePackingLimit::~bidispersePackingLimit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::kineticTheoryModels::packingLimitModels::
bidispersePackingLimit::alphaMax
(
    const label celli,
    const scalarList& ds
) const
{
    if (ds.size() != 2)
    {
        FatalErrorInFunction
            << typeName << " packing limit model only supports bidisperse "
            << "particle distributions, " << ds.size() << " in use."
            << exit(FatalError);
    }

    const phaseModel& phase1 = kt_.fluid().phases()[0];
    const phaseModel& phase2 = kt_.fluid().phases()[1];

    scalar alphaMax1 = phase1.alphaMax();
    scalar alpha1 = phase1[celli];
    scalar d1 = ds[0];

    scalar alphaMax2 = phase2.alphaMax();
    scalar alpha2 = phase2[celli];
    scalar d2 = ds[1];

    scalar Xi =
        max
        (
            alpha1/max(alpha1 + alpha2, residualAlpha_),
            alphaMax1/(alphaMax1 + (1.0 - alphaMax1)*alphaMax2)
        );

    scalar dByd = d1/d2;
    if (d1 > d2)
    {
        dByd = d2/d1;
    }

    return
        alphaMax1
      + (1.0 - sqrt(dByd))
       *(
            alphaMax1
          + (1.0 - alphaMax1)*alphaMax2
        )*(1.0 - Xi);
}


bool Foam::kineticTheoryModels::packingLimitModels::
bidispersePackingLimit::read()
{
    return true;
}


// ************************************************************************* //
