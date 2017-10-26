/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "ChaoConductivity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace conductivityModels
{
    defineTypeNameAndDebug(Chao, 0);

    addToRunTimeSelectionTable
    (
        conductivityModel,
        Chao,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::Chao::Chao
(
    const dictionary& dict,
    const polydisperseKineticTheoryModel& kt
)
:
    conductivityModel(dict, kt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::Chao::
~Chao()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::conductivityModels::Chao::kappa
(
    const phaseModel& phase,
    const volScalarField& Theta,
    const volScalarField& g0,
    const volScalarField& rho1,
    const volScalarField& da,
    const dimensionedScalar& e
) const
{
    volScalarField kCoeff
    (
        IOobject
        (
            "kCoeff",
            Theta.time().timeName(),
            Theta.mesh()
        ),
        Theta.mesh(),
        dimensionedScalar("0", dimLength, 0.0)
    );

    forAll(kt_.phases(), phasei)
    {
        const word& name2(kt_.phases()[phasei]);
        const phaseModel& phase2 = kt_.fluid().phases()[name2];
        const scalar& eij(kt_.es()[phasePairKey(phase.name(), name2)]);

        volScalarField m0
        (
            Foam::constant::mathematical::pi/6.0
           *(phase.rho()*pow3(da) + phase2.rho()*pow3(phase2.d()))
        );
        volScalarField dij(0.5*(da + phase2.d()));

        kCoeff +=
            phase2*phase2.rho()/m0*(1.0 + eij)*kt_.gs0(phase, phase2)*pow4(dij);
    }

    return
        constant::mathematical::pi*2.0/5.0*rho1*sqrt(Theta)*kCoeff;

}


bool Foam::kineticTheoryModels::conductivityModels::Chao::read()
{
    return true;
}


// ************************************************************************* //
