/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 OpenFOAM Foundation
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

#include "ReactingSolidPhaseModel.H"
#include "phaseSystem.H"
#include "fvMatrix.H"
#include "localEulerDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel, class SolidReactionType, class GasReactionType>
Foam::ReactingSolidPhaseModel
<
    BasePhaseModel,
    SolidReactionType,
    GasReactionType
>::ReactingSolidPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index),
    chemistryPtr_(basicSolidChemistryModel::New(this->thermo_())),
    gasPhaseName_(fluid.subDict(phaseName).lookup("gasPhase")),
    gasThermo_
    (
        fluid.mesh().template lookupObject<GasReactionType>
        (
            this->thermo().phasePropertyName
            (
                basicThermo::dictName,
                gasPhaseName_
            )
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel, class SolidReactionType, class GasReactionType>
Foam::ReactingSolidPhaseModel
<
    BasePhaseModel,
    SolidReactionType,
    GasReactionType
>::~ReactingSolidPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel, class SolidReactionType, class GasReactionType>
void Foam::ReactingSolidPhaseModel
<
    BasePhaseModel,
    SolidReactionType,
    GasReactionType
>::correctThermo()
{
    BasePhaseModel::correctThermo();
    chemistryPtr_->solve(this->fluid().mesh().time().deltaTValue());
}


template<class BasePhaseModel, class SolidReactionType, class GasReactionType>
Foam::tmp<Foam::fvScalarMatrix>
Foam::ReactingSolidPhaseModel
<
    BasePhaseModel,
    SolidReactionType,
    GasReactionType
>::R
(
    volScalarField& Yi
) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Yi, dimMass/dimTime));
    fvScalarMatrix& Su = tSu.ref();
    volScalarField RR
    (
        IOobject
        (
            IOobject::groupName("RR", Yi.name()),
            this->fluid().mesh().time().timeName(),
            this->fluid().mesh()
        ),
        this->fluid().mesh(),
        dimensionedScalar("0", dimDensity/dimTime, 0.0)
    );

    if
    (
        refCast<const SolidReactionType>
        (
            this->thermo()
        ).composition().contains(Yi.member())
    )
    {
        const label specieI =
            refCast<const SolidReactionType>
            (
                this->thermo()
            ).composition().species()[Yi.member()];

            RR.ref() = chemistryPtr_->RRs(specieI);
    }
    else if (chemistryPtr_->gasTable().contains(Yi.member()))
    {
        const label specieI =
            chemistryPtr_->gasTable()[Yi.member()];
        RR.ref() += chemistryPtr_->RRg(specieI)*(*this);
    }

    Su += RR;
    return tSu;
}


template<class BasePhaseModel, class SolidReactionType, class GasReactionType>
Foam::tmp<Foam::volScalarField>
Foam::ReactingSolidPhaseModel
<
    BasePhaseModel,
    SolidReactionType,
    GasReactionType
>::Qdot(const bool local) const
{
    tmp<volScalarField> tQdot
    (
        new volScalarField
        (
            IOobject
            (
                this->thermo().phasePropertyName("solidChemistry:Qdot"),
                this->fluid().mesh().time().timeName(),
                this->fluid().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->fluid().mesh(),
            dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );

    tQdot.ref() = chemistryPtr_->Qdot();
    if (!local)
    {
        tQdot.ref() *= -(*this);
    }

    return tQdot;
}


template<class BasePhaseModel, class SolidReactionType, class GasReactionType>
Foam::tmp<Foam::volScalarField>
Foam::ReactingSolidPhaseModel
<
    BasePhaseModel,
    SolidReactionType,
    GasReactionType
>::dmdt(const bool local) const
{
    tmp<volScalarField> tdmdt
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName
                (
                    "dmdt",
                    (local == true ? this->name() : gasPhaseName_)
                ),
                this->fluid().mesh().time().timeName(),
                this->fluid().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->fluid().mesh(),
            dimensionedScalar("dmdt", dimDensity/dimTime, 0.0)
        )
    );
    volScalarField& dmdt = tdmdt.ref();

    if (local)
    {
        label nSpecies =
            this->thermo_->composition().species().size();
        for (label specieI = 0; specieI < nSpecies; specieI++)
        {
            dmdt.ref() += chemistryPtr_->RRs(specieI);
        }
    }
    else
    {
        forAll(chemistryPtr_->gasTable(), specieI)
        {
            dmdt.ref() += chemistryPtr_->RRg(specieI);
        }
    }
    dmdt *= (*this);
    return tdmdt;
}

// ************************************************************************* //
