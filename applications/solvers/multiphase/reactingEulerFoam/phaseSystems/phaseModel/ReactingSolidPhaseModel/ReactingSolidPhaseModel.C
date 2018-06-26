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

template<class BasePhaseModel, class ReactionType>
Foam::ReactingSolidPhaseModel<BasePhaseModel, ReactionType>::
ReactingSolidPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index),
    chemistryPtr_(BasicChemistryModel<ReactionType>::New(this->thermo_())),
    active_(chemistryPtr_->lookupOrDefault("active", true)),
    integrateReactionRate_
    (
        chemistryPtr_->lookupOrDefault("integrateReactionRate", true)
    )
{
    if (integrateReactionRate_)
    {
        Info<< "    using integrated reaction rate" << endl;
    }
    else
    {
        Info<< "    using instantaneous reaction rate" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel, class ReactionType>
Foam::ReactingSolidPhaseModel<BasePhaseModel, ReactionType>::
~ReactingSolidPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel, class ReactionType>
void Foam::ReactingSolidPhaseModel<BasePhaseModel, ReactionType>::correctThermo()
{
    BasePhaseModel::correctThermo();

    if (active_)
    {
        if (integrateReactionRate_)
        {
            if (fv::localEulerDdt::enabled(this->fluid().mesh()))
            {
                const scalarField& rDeltaT =
                    fv::localEulerDdt::localRDeltaT(this->fluid().mesh());

                if (chemistryPtr_->found("maxIntegrationTime"))
                {
                    scalar maxIntegrationTime
                    (
                        readScalar
                        (
                            chemistryPtr_->lookup("maxIntegrationTime")
                        )
                    );

                    chemistryPtr_->solve
                    (
                        min(1.0/rDeltaT, maxIntegrationTime)()
                    );
                }
                else
                {
                    chemistryPtr_->solve((1.0/rDeltaT)());
                }
            }
            else
            {
                chemistryPtr_->solve
                (
                    this->fluid().mesh().time().deltaTValue()
                );
            }
        }
        else
        {
            chemistryPtr_->calculate();
        }
    }
}


template<class BasePhaseModel, class ReactionType>
Foam::tmp<Foam::fvScalarMatrix>
Foam::ReactingSolidPhaseModel<BasePhaseModel, ReactionType>::R
(
    volScalarField& Yi
) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Yi, dimMass/dimTime));

    fvScalarMatrix& Su = tSu.ref();

    if (active_)
    {
        const label specieI =
            refCast<const ReactionType>
            (
                this->thermo()
            ).composition().species()[Yi.member()];

        Su += chemistryPtr_->RR(specieI);
    }

    return tSu;
}


template<class BasePhaseModel, class ReactionType>
Foam::tmp<Foam::volScalarField>
Foam::ReactingSolidPhaseModel<BasePhaseModel, ReactionType>::Qdot() const
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

    if (active_)
    {
        tQdot.ref() = chemistryPtr_->Qdot();
    }

    return tQdot;
}


// ************************************************************************* //
