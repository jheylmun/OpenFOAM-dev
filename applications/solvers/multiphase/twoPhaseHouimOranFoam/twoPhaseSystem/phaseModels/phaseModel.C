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

#include "phaseModel.H"
#include "twoPhaseSystem.H"
#include "diameterModel.H"
#include "fvMatrix.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "dragModel.H"
#include "heatTransferModel.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "fixedValueFvsPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseModel, 0);
    defineRunTimeSelectionTable(phaseModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("alpha", dimless, 0)
    ),
    fluid_(fluid),
    name_(phaseName),
    phaseDict_
    (
        phaseProperties.subDict(name_)
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        fluid.subDict(phaseName).lookup("residualAlpha")
    ),
    residualRho_
    (
        "residualRho",
        dimDensity,
        fluid.subDict(phaseName).lookup("residualRho")
    ),
    thermoPtr_
    (
        rhoThermo::New(fluid.mesh(), phaseName)
    ),
    rho_(thermoPtr_->rho()),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),
    he_(thermoPtr_->he()),
    E_
    (
        IOobject
        (
            IOobject::groupName("E", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        he_ + 0.5*magSqr(U_),
        thermoPtr_->T().boundaryField().types()
    ),
    alphaRho_
    (
        IOobject
        (
            IOobject::groupName("alphaRho", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        (*this)*rho(),
        this->boundaryField().types()
    ),
    alphaRhoU_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoU", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        alphaRho_*U_,
        U_.boundaryField().types()
    ),
    alphaRhoE_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoE", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        alphaRho_*E_,
        E_.boundaryField().types()
    ),
    massFlux_
    (
        IOobject
        (
            IOobject::groupName("massFlux", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fvc::flux(U_)*fvc::interpolate((*this)*rho())
    ),
    momentumFlux_
    (
        IOobject
        (
            IOobject::groupName("momentumFlux", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        massFlux_*fvc::interpolate(U_)
    ),
    energyFlux_
    (
        IOobject
        (
            IOobject::groupName("energyFlux", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        massFlux_*fvc::interpolate(E_)
    ),
    gradAlpha_
    (
        IOobject
        (
            IOobject::groupName("gradAlpha", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fvc::grad(*this)
    ),
    fluxFunction_
    (
        phaseFluxFunction::New(fluid.mesh(), phaseName)
    )
{
    thermoPtr_->validate("phaseModel " + name_, "h", "e");

    const word phiName = IOobject::groupName("phi", name_);

    IOobject phiHeader
    (
        phiName,
        fluid_.mesh().time().timeName(),
        fluid_.mesh(),
        IOobject::NO_READ
    );

    dPtr_ = diameterModel::New
    (
        phaseDict_,
        *this
    );

    turbulence_ =
        PhaseCompressibleTurbulenceModel<phaseModel>::New
        (
            *this,
            this->rho_,
            this->U_,
            this->massFlux_,
            phi(),
            *this
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseModel::~phaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseModel& Foam::phaseModel::otherPhase() const
{
    return fluid_.otherPhase(*this);
}


const Foam::PhaseCompressibleTurbulenceModel<Foam::phaseModel>&
Foam::phaseModel::turbulence() const
{
    return turbulence_();
}

Foam::PhaseCompressibleTurbulenceModel<Foam::phaseModel>&
Foam::phaseModel::turbulence()
{
    return turbulence_();
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::d() const
{
    return dPtr_().d();
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::nuEff() const
{
    return turbulence_->nuEff();
}


Foam::tmp<Foam::scalarField> Foam::phaseModel::nuEff(const label patchi) const
{
    return turbulence_->nuEff(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::k() const
{
    return turbulence_->k();
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::epsilon() const
{
    return turbulence_->epsilon();
}


Foam::tmp<Foam::volSymmTensorField> Foam::phaseModel::R() const
{
    return turbulence_->R();
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::pPrime() const
{
    return turbulence_->pPrime();
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseModel::pPrimef() const
{
    return turbulence_->pPrimef();
}


Foam::tmp<Foam::volSymmTensorField> Foam::phaseModel::devRhoReff() const
{
    return turbulence_->devRhoReff();
}


Foam::tmp<Foam::fvVectorMatrix> Foam::phaseModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return turbulence_->divDevRhoReff(U);
}


void Foam::phaseModel::correct()
{
    return dPtr_->correct();
}


void Foam::phaseModel::setNSteps(const label nSteps)
{
    bool setAlpha = !otherPhase().granular();
    if (setAlpha)
    {
        alphas_.setSize(nSteps);
        deltaAlpha_.setSize(nSteps);
    }
    alphaRhos_.setSize(nSteps);
    deltaAlphaRho_.setSize(nSteps);
    alphaRhoUs_.setSize(nSteps);
    deltaAlphaRhoU_.setSize(nSteps);
    alphaRhoEs_.setSize(nSteps);
    deltaAlphaRhoE_.setSize(nSteps);

    if (granular())
    {
        alphaRhoPTEs_.setSize(nSteps);
        deltaAlphaRhoPTE_.setSize(nSteps);
    }

    forAll (deltaAlphaRho_, stepi)
    {
        if (setAlpha)
        {
            alphas_.set
            (
                stepi,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "alpha"+Foam::name(stepi),
                            name()
                        ),
                        fluid_.mesh().time().timeName(),
                        fluid_.mesh()
                    ),
                    fluid_.mesh(),
                    dimensionedScalar("zero", dimless, 0.0)
                )
            );
            deltaAlpha_.set
            (
                stepi,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "deltaAlpha"+Foam::name(stepi),
                            name()
                        ),
                        fluid_.mesh().time().timeName(),
                        fluid_.mesh()
                    ),
                    fluid_.mesh(),
                    dimensionedScalar("zero", dimless/dimTime, 0.0)
                )
            );
        }
        alphaRhos_.set
        (
            stepi,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "alphaRho" + Foam::name(stepi),
                        name()
                    ),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh()
                ),
                fluid_.mesh(),
                dimensionedScalar("zero", alphaRho_.dimensions(), 0.0)
            )
        );
        deltaAlphaRho_.set
        (
            stepi,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "deltaAlphaRho" + Foam::name(stepi),
                        name()
                    ),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh()
                ),
                fluid_.mesh(),
                dimensionedScalar("zero", alphaRho_.dimensions()/dimTime, 0.0)
            )
        );
        alphaRhoUs_.set
        (
            stepi,
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "alphaRhoU" + Foam::name(stepi),
                        name()
                    ),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh()
                ),
                fluid_.mesh(),
                dimensionedVector("zero", alphaRhoU_.dimensions(), Zero)
            )
        );
        deltaAlphaRhoU_.set
        (
            stepi,
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "deltaAlphaRhoU" + Foam::name(stepi),
                        name()
                    ),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh()
                ),
                fluid_.mesh(),
                dimensionedVector("zero", alphaRhoU_.dimensions()/dimTime, Zero)
            )
        );
        alphaRhoEs_.set
        (
            stepi,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "alphaRhoE" + Foam::name(stepi),
                        name()
                    ),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh()
                ),
                fluid_.mesh(),
                dimensionedScalar("zero", alphaRhoE_.dimensions(), 0.0)
            )
        );
        deltaAlphaRhoE_.set
        (
            stepi,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "deltaAlphaRhoE" + Foam::name(stepi),
                        name()
                    ),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh()
                ),
                fluid_.mesh(),
                dimensionedScalar("zero", alphaRhoE_.dimensions()/dimTime, 0.0)
            )
        );

        if (granular())
        {
            alphaRhoPTEs_.set
            (
                stepi,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "alphaRhoPTE" + Foam::name(stepi),
                            name()
                        ),
                        fluid_.mesh().time().timeName(),
                        fluid_.mesh()
                    ),
                    fluid_.mesh(),
                    dimensionedScalar("zero", alphaRhoE_.dimensions(), 0.0)
                )
            );
            deltaAlphaRhoPTE_.set
            (
                stepi,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "deltaAlphaRhoPTE" + Foam::name(stepi),
                            name()
                        ),
                        fluid_.mesh().time().timeName(),
                        fluid_.mesh()
                    ),
                    fluid_.mesh(),
                    dimensionedScalar
                    (
                        "zero",
                        alphaRhoE_.dimensions()/dimTime,
                        0.0
                    )
                )
            );
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::dissipationCoeff()const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ProductionCoeff",
                fluid_.mesh().time().timeName(),
                fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fluid_.mesh(),
            dimensionedScalar("zero", dimensionSet(1, -2, -2, 0, 0), 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::productionCoeff()const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "productionCoeff",
                fluid_.mesh().time().timeName(),
                fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fluid_.mesh(),
            dimensionedScalar("zero", dimensionSet(1, -2, -2, 0, 0), 0.0)
        )
    );
}


bool Foam::phaseModel::read(const dictionary& phaseProperties)
{
    phaseDict_ = phaseProperties.subDict(name_);
    return dPtr_->read(phaseDict_);
}


void Foam::phaseModel::correctInflowOutflow(surfaceScalarField& alphaPhi) const
{
    surfaceScalarField::Boundary& alphaPhiBf = alphaPhi.boundaryFieldRef();
    const volScalarField::Boundary& alphaBf = boundaryField();
    const surfaceScalarField::Boundary& phiBf = phi().boundaryField();

    forAll(alphaPhiBf, patchi)
    {
        fvsPatchScalarField& alphaPhip = alphaPhiBf[patchi];

        if (!alphaPhip.coupled())
        {
            alphaPhip = phiBf[patchi]*alphaBf[patchi];
        }
    }
}


// ************************************************************************* //
