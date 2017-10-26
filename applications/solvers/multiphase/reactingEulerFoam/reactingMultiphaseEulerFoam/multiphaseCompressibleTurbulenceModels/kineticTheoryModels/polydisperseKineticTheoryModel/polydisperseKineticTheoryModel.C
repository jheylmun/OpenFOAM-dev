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
#include "multiphaseSystem.H"
#include "radialModel.H"
#include "viscosityModel.H"
#include "conductivityModel.H"
#include "granularPressureModel.H"
#include "frictionalStressModel.H"
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
    phases_(),
    radialModel_
    (
        kineticTheoryModels::radialModel::New
        (
            dict_,
            *this
        )
    ),
    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            dict_,
            *this
        )
    ),
    conductivityModel_
    (
        kineticTheoryModels::conductivityModel::New
        (
            dict_,
            *this
        )
    ),
    granularPressureModel_
    (
        kineticTheoryModels::granularPressureModel::New
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
    PsCoeffs_(),
    PsSrcs_(),
    eTable_
    (
        dict_.lookup("coeffRest")
    ),
    CfTable_
    (
        dict_.lookup("coeffFric")
    ),
    alphap_
    (
        IOobject
        (
            "alphap",
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimless, 0.0)
    ),
    Up_
    (
        IOobject
        (
            "Up",
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedVector("0", dimVelocity, Zero)
    ),
    alphaMax_
    (
        IOobject
        (
            "packingLimit",
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
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
        granularPressureModel_->read();
        viscosityModel_->read();
        conductivityModel_->read();
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


Foam::tmp<Foam::volScalarField> Foam::polydisperseKineticTheoryModel::nu
(
    const phaseModel& phase
) const
{
    return viscosityModel_->nu
    (
        phase,
        ktModels_[phase.name()]->Theta(),
        ktModels_[phase.name()]->gs0(),
        phase.rho(),
        phase.d(),
        eTable_[phasePairKey(phase.name(), phase.name())]
    );
}


Foam::tmp<Foam::volScalarField> Foam::polydisperseKineticTheoryModel::kappa
(
    const phaseModel& phase
) const
{
    return conductivityModel_->kappa
    (
        phase,
        ktModels_[phase.name()]->Theta(),
        ktModels_[phase.name()]->gs0(),
        phase.rho(),
        phase.d(),
        eTable_[phasePairKey(phase.name(), phase.name())]
    );
}


Foam::tmp<Foam::volScalarField>
Foam::polydisperseKineticTheoryModel::
granularPressureCoeff(const phaseModel& phase) const
{
    return *PsCoeffs_[phase.name()];
}


Foam::tmp<Foam::volScalarField>
Foam::polydisperseKineticTheoryModel::
granularPressureSrc(const phaseModel& phase) const
{
    return *PsSrcs_[phase.name()];
}


Foam::tmp<Foam::volScalarField>
Foam::polydisperseKineticTheoryModel::
granularPressurePrimeCoeff(const phaseModel& phase) const
{
    return *PsPrimeCoeffs_[phase.name()];
}


Foam::tmp<Foam::volScalarField>
Foam::polydisperseKineticTheoryModel::
granularPressurePrimeSrc(const phaseModel& phase) const
{
    return *PsPrimeSrcs_[phase.name()];
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
    ktModels_.set
    (
        phaseName,
        const_cast<RASModels::kineticTheoryModel*>(&ktModel)
    );

    PsCoeffs_.set
    (
        phaseName,
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("PsCoeff", phaseName),
                fluid_.mesh().time().timeName(),
                fluid_.mesh()
            ),
            fluid_.mesh(),
            dimensionedScalar("PsCoeff", dimDensity, 0.0)
        )
    );
    PsPrimeCoeffs_.set
    (
        phaseName,
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("PsPrimeCoeff", phaseName),
                fluid_.mesh().time().timeName(),
                fluid_.mesh()
            ),
            fluid_.mesh(),
            dimensionedScalar("PsPrimeCoeff", dimDensity, 0.0)
        )
    );
    PsSrcs_.set
    (
        phaseName,
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("PsSrc", phaseName),
                fluid_.mesh().time().timeName(),
                fluid_.mesh()
            ),
            fluid_.mesh(),
            dimensionedScalar("PsSrc", dimPressure, 0.0)
        )
    );
    PsPrimeSrcs_.set
    (
        phaseName,
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("PsPrimeSrc", phaseName),
                fluid_.mesh().time().timeName(),
                fluid_.mesh()
            ),
            fluid_.mesh(),
            dimensionedScalar("PsPrimeSrc", dimPressure, 0.0)
        )
    );
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
    Up_ /= max(alphap_, 1e-6);

    // Set packing limit
    forAll(alphap_, celli)
    {
        SortableList<scalar> ds(phases_.size());
        forAll(phases_, phasei)
        {
            ds[phasei] = fluid_.phases()[phases_[phasei]].d()()[celli];
        }

        ds.reverseSort();
        alphaMax_[celli] = 0.63;//calcAlphaMax(0, celli, ds);
    }
    alphaMax_.correctBoundaryConditions();

    forAll(phases_, phasei)
    {
        const word& name(phases_[phasei]);
        *PsCoeffs_[name] = dimensionedScalar("0", dimDensity, 0.0);
        *PsPrimeCoeffs_[name] = dimensionedScalar("0", dimDensity, 0.0);
        *PsSrcs_[name] =
            dimensionedScalar("0", dimPressure, 0.0);
        *PsPrimeSrcs_[name] =
            dimensionedScalar("0", dimPressure, 0.0);

        const phaseModel& phase(fluid_.phases()[name]);
        const volScalarField& rho1(phase.rho());
        volScalarField d1(phase.d());

        forAll(phases_, phasej)
        {
            const word& name2(phases_[phasej]);
            const phaseModel& phase2(fluid_.phases()[name2]);
            const volScalarField& rho2(phase2.rho());
            const volScalarField& Theta2(ktModels_[name2]->Theta());
            volScalarField d2(phase2.d());

            const scalar& eij(eTable_[phasePairKey(name, name2, false)]);
            tmp<volScalarField> gij(gs0(phase, phase2));
            tmp<volScalarField> gijPrime(gs0Prime(phase, phase2));

            volScalarField m0
            (
                constant::mathematical::pi/6.0
               *(rho1*pow3(d1) + rho2*pow3(d2))
            );

            volScalarField PsCoeff
            (
                granularPressureModel_->granularPressureCoeff
                (
                    phase,
                    phase2,
                    gij,
                    eij
                )
            );
            volScalarField PsPrimeCoeff
            (
                granularPressureModel_->granularPressureCoeffPrime
                (
                    phase,
                    phase2,
                    gij,
                    gijPrime,
                    eij
                )
            );

            *PsPrimeCoeffs_[name] += PsPrimeCoeff;
            *PsCoeffs_[name] += PsCoeff;

            *PsPrimeSrcs_[name] +=
                PsPrimeCoeff*(Theta2 + 0.2*magSqr(phase.U() - phase2.U()));

            *PsSrcs_[name] +=
                PsCoeff*(Theta2 + 0.2*magSqr(phase.U() - phase2.U()));
        }
    }
}


// ************************************************************************* //
