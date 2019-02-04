/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "dragOde.H"

void Foam::dragOde::setUs()
{
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase = phaseModels_[phasei];
        if (!phase.stationary())
        {
            label ui = phase.index()*nDims_;
            if (solutionD_[0] == 1)
            {
                phase.URef()[celli_].x() = q_[ui];
                ui++;
            }
            if (solutionD_[1] == 1)
            {
                phase.URef()[celli_].y() = q_[ui];
                ui++;

            }
            if (solutionD_[2] == 1)
            {
                phase.URef()[celli_].z() = q_[ui];
            }
        }
    }
}


void Foam::dragOde::setq()
{
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase = phaseModels_[phasei];
        label ui = phase.index()*nDims_;
        if (solutionD_[0] == 1)
        {
            q_[ui] = phase.URef()[celli_].x();
            ui++;
        }
        if (solutionD_[1] == 1)
        {
            q_[ui] = phase.URef()[celli_].y();
            ui++;
        }
        if (solutionD_[2] == 1)
        {
            q_[ui] = phase.URef()[celli_].z();
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragOde::dragOde(phaseSystem& fluid, dragModelTable& dragModels)
:
    ODESystem(),
    dict_(fluid.subDict("dragOdeCoeffs")),
    solveDrag_(dict_.lookupOrDefault("solveOde", false)),
    fluid_(fluid),
    phaseModels_(fluid.phases()),
    solutionD_(fluid.mesh().solutionD()),
    nDims_(fluid_.mesh().nSolutionD()),
    nEqns_(nDims_*phaseModels_.size()),
    q_(nEqns_, 0.0),
    dqdt_(nEqns_, 0.0),
    deltaT_(phaseModels_[0].size(), fluid.mesh().time().deltaTValue())
{
    if (solveDrag_)
    {
        odeSolver_ = ODESolver::New(*this, dict_);
    }

    //- Add unorded phase pairs with vaild drag models
    forAllIter
    (
        dragModelTable,
        dragModels,
        dragModelIter
    )
    {
        const phasePair& pair(fluid.phasePairs()[dragModelIter.key()]);

        const phasePairKey key(pair.first(), pair.second());

        phasePairs_.append(new phasePair(pair));
        dragModels_.append
        (
            &dragModelIter()()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragOde::~dragOde()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::dragOde::derivatives
(
    const scalar time,
    const scalarField& q,
    scalarField& dqdt
) const
{
    dqdt = 0.0;

    forAll(phasePairs_, pairi)
    {
        const phasePair& pair(phasePairs_[pairi]);
        scalar drag = dragModels_[pairi].cellK(celli_);

        const phaseModel& phase1(pair.phase1());
        const phaseModel& phase2(pair.phase2());
        label ui = phase1.index()*nDims_;
        scalar alphaRho1 =
            Foam::max
            (
                phase1[celli_],
                phase1.residualAlpha().value()
            )*phase1.thermo().cellrho(celli_);
        scalar u1 = q_[ui];

        label uj = phase2.index()*nDims_;
        scalar alphaRho2 =
            Foam::max
            (
                phase2[celli_],
                phase2.residualAlpha().value()
            )*phase2.thermo().cellrho(celli_);
        scalar u2 = q_[uj];

        scalar drag1 = drag/alphaRho1;
        scalar drag2 = drag/alphaRho2;

        if (nDims_ > 0)
        {
            dqdt[ui] += drag1*(u2 - u1);
            dqdt[uj] += drag2*(u1 - u2);
            ui++;
            uj++;
        }
        if (nDims_ > 1)
        {
            dqdt[ui] += drag1*(u2 - u1);
            dqdt[uj] += drag2*(u1 - u2);
            ui++;
            uj++;
        }
        if (nDims_ > 2)
        {
            dqdt[ui] += drag1*(u2 - u1);
            dqdt[uj] += drag2*(u1 - u2);
        }
    }
}


void Foam::dragOde::jacobian
(
    const scalar t,
    const scalarField& q,
    scalarField& dqdt,
    scalarSquareMatrix& J
) const
{
    dqdt = 0.0;
    J = scalarSquareMatrix(nEqns_, 0.0);
    forAll(phasePairs_, pairi)
    {
        const phasePair& pair(phasePairs_[pairi]);
        scalar drag = dragModels_[pairi].cellK(celli_);

        const phaseModel& phase1(pair.phase1());
        const phaseModel& phase2(pair.phase2());
        label ui = phase1.index()*nDims_;
        scalar alphaRho1 =
            Foam::max
            (
                phase1[celli_],
                phase1.residualAlpha().value()
            )*phase1.thermo().cellrho(celli_);
        scalar u1 = q_[ui];

        label uj = phase2.index()*nDims_;
        scalar alphaRho2 =
            Foam::max
            (
                phase2[celli_],
                phase2.residualAlpha().value()
            )*phase2.thermo().cellrho(celli_);
        scalar u2 = q_[uj];

        scalar drag1 = drag/alphaRho1;
        scalar drag2 = drag/alphaRho2;

        if (nDims_ > 0)
        {
            dqdt[ui] += drag1*(u2 - u1);
            dqdt[uj] += drag2*(u1 - u2);

            J(ui, ui) += -drag1;
            J(ui, uj) += drag1;
            J(uj, uj) += -drag2;
            J(uj, ui) += drag2;

            ui++;
            uj++;
        }
        if (nDims_ > 1)
        {
            dqdt[ui] += drag1*(u2 - u1);
            dqdt[uj] += drag2*(u1 - u2);

            J(ui, ui) += -drag1;
            J(ui, uj) += drag1;
            J(uj, uj) += -drag2;
            J(uj, ui) += drag2;

            ui++;
            uj++;
        }
        if (nDims_ > 2)
        {
            dqdt[ui] += drag1*(u2 - u1);
            dqdt[uj] += drag2*(u1 - u2);
            J(ui, ui) += -drag1;
            J(ui, uj) += drag1;
            J(uj, uj) += -drag2;
            J(uj, ui) += drag2;
        }
    }
}


Foam::scalar Foam::dragOde::solve
(
    const scalar& deltaT
)
{
    for(celli_ = 0; celli_ < phaseModels_[0].size(); celli_++)
    {
        setq();

        scalar timeLeft = deltaT;
        while (timeLeft > small)
        {
            scalar dt = timeLeft;
            odeSolver_->solve(0, dt, q_, deltaT_[celli_]);
            setUs();
            timeLeft -= dt;
            deltaT_[celli_] = dt;
        }
    }
    return min(deltaT_);
}




// ************************************************************************* //
