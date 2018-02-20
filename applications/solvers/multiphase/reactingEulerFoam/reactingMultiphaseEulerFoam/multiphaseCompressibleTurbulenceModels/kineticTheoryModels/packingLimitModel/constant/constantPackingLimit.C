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

#include "constantPackingLimit.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace packingLimitModels
{
    defineTypeNameAndDebug(constantPackingLimit, 0);

    addToRunTimeSelectionTable
    (
        packingLimitModel,
        constantPackingLimit,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::packingLimitModels::
constantPackingLimit::constantPackingLimit
(
    const dictionary& dict,
    const polydisperseKineticTheoryModel& kt
)
:
    packingLimitModel(dict, kt),
    maxAlpha_(dict.lookupType<scalar>("alphaMax"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::packingLimitModels::
constantPackingLimit::~constantPackingLimit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::kineticTheoryModels::packingLimitModels::
constantPackingLimit::alphaMax
(
    const label celli,
    const scalarList& ds
) const
{
    return maxAlpha_;
}


bool Foam::kineticTheoryModels::packingLimitModels::
constantPackingLimit::read()
{
    maxAlpha_ = dict_.lookupType<scalar>("alphaMax");

    return true;
}


// ************************************************************************* //
