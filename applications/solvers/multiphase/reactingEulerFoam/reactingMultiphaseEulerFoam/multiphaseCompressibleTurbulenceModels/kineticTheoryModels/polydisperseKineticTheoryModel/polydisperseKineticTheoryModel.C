/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "polydisperseKineticTheoryModel.H"
#include "kineticTheoryModel.H"
#include "multiphaseSystem.H"
#include "radialModel.H"
#include "frictionalStressModel.H"
// #include "granularPressureModel.H"
#include "phaseSystem.H"
#include "mathematicalConstants.H"
#include "SortableList.H"
#include "zeroGradientFvPatchFields.H"


// * * * * * * * * * * * * * * *  Private Funtions * * * * * * * * * * * * * //

Foam::scalar Foam::polydisperseKineticTheoryModel::calcAlphaMax
(
    const label phasei,
    const label cellI,
    const scalarList& ds
) const
{
    if (phasei + 1 == phases_.size())
    {
        return fluid_.phases()[phasei].alphaMax();
    }
    const phaseModel& phase1 = fluid_.phases()[phasei];
    const label phaseJ(phasei + 1);
    const phaseModel& phase2 = fluid_.phases()[phaseJ];

    scalar alphaMax2 = calcAlphaMax(phaseJ, cellI, ds);

    return
        phase1.alphaMax()
      + (1.0 - sqrt(ds[phaseJ]/ds[phasei]))
       *(
            phase1.alphaMax()
          + (1.0 - phase1.alphaMax())*alphaMax2
        )*phase2[cellI]/max(alphap_[cellI], 1e-6);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polydisperseKineticTheoryModel::polydisperseKineticTheoryModel
(
    const phaseSystem& fluid
)
:
    regIOobject
    (
        IOobject
        (
            "polydisperseKineticTheory",
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    fluid_(fluid),
    dict_(fluid.subDict("kineticTheory")),
    name_(dict_.lookup("name")),
    phases_(),
    radialModel_
    (
        kineticTheoryModels::radialModel::New
        (
            dict_,
            *this
        )
    ),
    frictionalStressModel_
    (
        kineticTheoryModels::frictionalStressModel::New
        (
            dict_,
            *this
        )
    ),
//     granularPressureModel_
//     (
//         kineticTheoryModels::granularPressureModel::New
//         (
//             dict_,
//             *this
//         )
//     ),
    eTable_
    (
        dict_.lookup("coeffRest")
    ),
    CfTable_
    (
        dict_.lookup("coeffFric")
    ),
    alphaMax_
    (
        IOobject
        (
            IOobject::groupName("alphaMax", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimless, 0.0),
        wordList
        (
            alphap_.boundaryField().size(),
            zeroGradientFvPatchScalarField::typeName
        )
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict_
    ),
    alphap_
    (
        IOobject
        (
            IOobject::groupName("alpha", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimless, 0.0)
    ),
    Up_
    (
        IOobject
        (
            IOobject::groupName("U", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedVector("0", dimVelocity, Zero)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polydisperseKineticTheoryModel::~polydisperseKineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::polydisperseKineticTheoryModel::read()
{
    if (fluid_.modified())
    {
        residualAlpha_.readIfPresent(dict_);

        radialModel_->read();
        frictionalStressModel_->read();

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField> Foam::polydisperseKineticTheoryModel::gs0
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return radialModel_->g0(phase1, phase2);
}


Foam::tmp<Foam::volScalarField> Foam::polydisperseKineticTheoryModel::gs0Prime
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return radialModel_->g0prime(phase1, phase2);
}


Foam::tmp<Foam::volScalarField>
Foam::polydisperseKineticTheoryModel::
frictionalPressure(const phaseModel& phase) const
{
    return frictionalStressModel_->frictionalPressure
    (
        phase,
        alphaMax_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::polydisperseKineticTheoryModel::
frictionalPressurePrime(const phaseModel& phase) const
{
    return frictionalStressModel_->frictionalPressurePrime
    (
        phase,
        alphaMax_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::polydisperseKineticTheoryModel::nuFrictional(const phaseModel& phase) const
{
    tmp<volTensorField> tgradU(fvc::grad(Up_));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));

    return frictionalStressModel_->nu
    (
        phase,
        alphaMax_,
        frictionalPressure(phase)/phase.rho(),
        D
    );
}



const Foam::wordList& Foam::polydisperseKineticTheoryModel::phases() const
{
    return phases_;
}


void Foam::polydisperseKineticTheoryModel::addPhase
(
    const RASModels::kineticTheoryModel& ktModel
)
{
    const word& phaseName(ktModel.phase().name());
    phases_.append(phaseName);

    forAll(phases_, phasei)
    {
        phasePairKey key
        (
            phaseName,
            phases_[phasei],
            false
        );
        pairs_.insert
        (
            key,
            autoPtr<phasePair>
            (
                new phasePair
                (
                    fluid_.phases()[phaseName],
                    fluid_.phases()[phases_[phasei]]
                )
            )
        );

        PsCoeffs_.insert
        (
            key,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("PsCoeff", pairs_[key]().name()),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh()
                ),
                fluid_.mesh(),
                dimensionedScalar("PsCoeff", dimless, 0.0)
            )
        );

    }
}


bool Foam::polydisperseKineticTheoryModel::found(const word& phaseName) const
{
    forAll(phases_, phasei)
    {
        if (phases_[phasei] == phaseName)
        {
            return true;
        }
    }
    return false;
}


void Foam::polydisperseKineticTheoryModel::correct()
{
    alphap_ = 0.0;
    Up_ = dimensionedVector("0", dimVelocity, Zero);

    forAll(phases_, phasei)
    {
        const volScalarField& alpha = fluid_.phases()[phases_[phasei]];
        alphap_ += alpha;
        Up_ += alpha*fluid_.phases()[phases_[phasei]].U();
    }
    Up_ /= max(alphap(), residualAlpha_);

    if (Switch(dict_.lookup("constantPackingLimit")))
    {
        alphaMax_ = dimensionedScalar("alphaMax", dimless, dict_);
    }
    else
    {
        // Set packing limit
        forAll(alphap_, celli)
        {
            SortableList<scalar> ds(phases_.size());
            forAll(phases_, phasei)
            {
                ds[phasei] = fluid_.phases()[phases_[phasei]].d()()[celli];
            }

            ds.reverseSort();
            alphaMax_[celli] = calcAlphaMax(0, celli, ds);
        }
    }
    alphaMax_.correctBoundaryConditions();
}


void Foam::polydisperseKineticTheoryModel::correctAlphap()
{
    alphap_ = 0;
    forAll(phases_, phasei)
    {
        alphap_ += fluid_.phases()[phases_[phasei]];
    }
}
// ************************************************************************* //
