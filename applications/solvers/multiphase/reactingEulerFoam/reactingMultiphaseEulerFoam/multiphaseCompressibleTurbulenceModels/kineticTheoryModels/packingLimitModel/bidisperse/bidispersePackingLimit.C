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
    const labelList& indices,
    const scalarList& ds
) const
{
    label index1 = 0;
    label index2 = indices.size() - 1;

    const phaseModel& phase1 = kt_.fluid().phases()[indices[index1]];
    const phaseModel& phase2 = kt_.fluid().phases()[indices[index2]];

    scalar alphaMax2 = phase2.alphaMax();

    return
        phase1.alphaMax()
      + (1.0 - sqrt(ds[index2]/ds[index1]))
       *(
            phase1.alphaMax()
          + (1.0 - phase1.alphaMax())*alphaMax2
        )*phase2[celli]/max(kt_.alphap()[celli], residualAlpha_);
}


bool Foam::kineticTheoryModels::packingLimitModels::
bidispersePackingLimit::read()
{
    return true;
}


// ************************************************************************* //
