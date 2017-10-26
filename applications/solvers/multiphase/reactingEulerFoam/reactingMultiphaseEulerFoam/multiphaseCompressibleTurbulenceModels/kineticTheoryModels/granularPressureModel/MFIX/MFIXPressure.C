/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "GidaspowPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace granularPressureModels
{
    defineTypeNameAndDebug(Gidaspow, 0);

    addToRunTimeSelectionTable
    (
        granularPressureModel,
        Gidaspow,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Gidaspow::Gidaspow
(
    const dictionary& dict
)
:
    granularPressureModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Gidaspow::~Gidaspow()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Gidaspow::
granularPressureCoeff
(
    const word& phaseName
) const
{

    const phaseModel& phase(kt_.fluid().phases()[phaseName]);
    const volScalarField& alpha(phase);
    const volScalarField& rho(phase.rho();
    const volScalarField& d(phase.d());

    volScalarField PsCoeffij
    (
        IOobject
        (
            "PsCoeffij",
            alpha.time().timeName(),
            alpha.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
        alpha.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    forAll(kt_.phases(), phaseJ)
    {
        const phaseModel& phase2(kt_.fluid().phases()[kt_.phases()[phaseJ]);
        tmp<volScalarField> d12 = 0.5*(d + phase2.d());
        scalar eij(kt_.coeffRest()[phasePair(phase1.name(), phase2.name())]);
        tmp<volScalarField> g0(kt_.gs0(phaseName, kt_.phases()[phaseJ]));

        PsCoeffij += 2.0*pow3(d12/d)*(1.0 + eij)*g0*phase2;
    }

    return rho1*alpha1*PsCoeffij;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Gidaspow::
granularPressureCoeffPrime
(
    const word& phaseName
) const
{
    return rho1*(1.0 + alpha1*(1.0 + e)*(4.0*g0 + 2.0*g0prime*alpha1));
}


// ************************************************************************* //
