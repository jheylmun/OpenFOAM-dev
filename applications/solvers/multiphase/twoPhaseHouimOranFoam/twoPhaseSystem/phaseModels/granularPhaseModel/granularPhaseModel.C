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

#include "granularPhaseModel.H"
#include "twoPhaseSystem.H"
#include "diameterModel.H"
#include "fvMatrix.H"
#include "dragModel.H"
#include "heatTransferModel.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        granularPhaseModel,
        dictionary,
        granular
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::granularPhaseModel::granularPhaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    phaseModel(fluid, phaseProperties, phaseName),
    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            this->phaseDict_
        )
    ),
    conductivityModel_
    (
        kineticTheoryModels::conductivityModel::New
        (
            this->phaseDict_
        )
    ),
    radialModel_
    (
        kineticTheoryModels::radialModel::New
        (
            this->phaseDict_
        )
    ),
    granularPressureModel_
    (
        kineticTheoryModels::granularPressureModel::New
        (
            this->phaseDict_
        )
    ),
    frictionalStressModel_
    (
        kineticTheoryModels::frictionalStressModel::New
        (
            this->phaseDict_
        )
    ),

    equilibrium_(this->phaseDict_.lookup("equilibrium")),
    e_("e", dimless, this->phaseDict_),
    alphaMax_("alphaMax", dimless, this->phaseDict_),
    alphaMinFriction_
    (
        "alphaMinFriction",
        dimless,
        this->phaseDict_
    ),

    maxNut_
    (
        "maxNut",
        dimensionSet(0,2,-1,0,0),
        this->phaseDict_.lookupOrDefault<scalar>("maxNut",1000)
    ),

    Theta_
    (
        IOobject
        (
            IOobject::groupName("Theta", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),

    Ps_
    (
        IOobject
        (
            IOobject::groupName("Ps", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (
            Theta_
           *granularPressureModel_->granularPressureCoeff
            (
                (*this),
                radialModel_->g0(*this, alphaMinFriction_, alphaMax_),
                this->rho_,
                e_
            )
        )
    ),
    Pfric_
    (
        IOobject
        (
            IOobject::groupName("Pfric", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        frictionalStressModel_->frictionalPressure
        (
            *this,
            alphaMinFriction_,
            alphaMax_
        )
    ),

    lambda_
    (
        IOobject
        (
            IOobject::groupName("lambda", phaseName),
            fluid.time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),

    gs0_
    (
        IOobject
        (
            IOobject::groupName("gs0", phaseName),
            fluid.time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),

    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", phaseName),
            fluid.time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),

    nut_
    (
        IOobject
        (
            IOobject::groupName("nut", phaseName),
            fluid.time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),

    nuFric_
    (
        IOobject
        (
            IOobject::groupName("nuFric", phaseName),
            fluid.time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),
    alphaRhoPTE_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPTE", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        1.5*alphaRho_*Theta_,
        Theta_.boundaryField().types()
    ),
    PTEFlux_
    (
        IOobject
        (
            IOobject::groupName("PTEFlux", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        1.5*massFlux_*fvc::interpolate(Theta_)
    )
{
    // Kinetic energy is not included in E
    E_ = he_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularPhaseModel::~granularPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::granularPhaseModel::c() const
{
    volScalarField A(1.0 + 2.0*(1.0 + e_)*(*this)*gs0_);
    tmp<volScalarField> B
    (
        2.0*(1.0 + e_)
       *(
            gs0_
          + radialModel_->g0prime
            (
                *this,
                alphaMinFriction_,
                alphaMax_
            )*(*this)
        )
    );

    volScalarField dAlphaCrit(*this - alphaMinFriction_);
    volScalarField dAlphaMax(alphaMax_ - *this);

    tmp<volScalarField> cFric
    (
        (Ps_ + Pfric_)*dAlphaCrit/(pow(dAlphaMax, 5)*(this->rho_))
       *(
            (*this)*(0.2 + 0.5*dAlphaCrit/dAlphaMax)
          + 0.1*dAlphaCrit
        )
    );
    cFric.ref() *= pos(dAlphaCrit);

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("a", name()),
            sqrt(Theta_*(A + 2.0/3.0*sqr(A) + (*this)*B) + cFric)
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::granularPhaseModel::k() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField> Foam::granularPhaseModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField> Foam::granularPhaseModel::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("R", name()),
                this->fluid_.mesh().time().timeName(),
                this->fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - (nut_)*dev(twoSymm(fvc::grad(U())))
          - (lambda_*fvc::div(U()))*symmTensor::I
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::granularPhaseModel::pPrime() const
{
    tmp<volScalarField> tpPrime
    (
        Theta_
       *granularPressureModel_->granularPressureCoeffPrime
        (
            *this,
            radialModel_->g0(*this, alphaMinFriction_, alphaMax_),
            radialModel_->g0prime(*this, alphaMinFriction_, alphaMax_),
            this->rho_,
            e_
        )
     +  frictionalStressModel_->frictionalPressurePrime
        (
            *this,
            alphaMinFriction_,
            alphaMax_
        )
    );

    volScalarField::Boundary& bpPrime =
        tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi] == 0;
        }
    }

    return tpPrime;
}


Foam::tmp<Foam::surfaceScalarField> Foam::granularPhaseModel::pPrimef() const
{
    return fvc::interpolate(pPrime());
}


Foam::tmp<Foam::volSymmTensorField>
Foam::granularPhaseModel::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", this->name()),
                this->fluid_.mesh().time().timeName(),
                this->fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - (this->rho_*nut_)*dev(twoSymm(fvc::grad(this->U_)))
          - ((this->rho_*lambda_)*fvc::div(this->phiPtr_()))*symmTensor::I
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix> Foam::granularPhaseModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
      - fvm::laplacian(rho()*nut_, U)
      - fvc::div
        (
            (this->rho_*nut_)*dev2(::Foam::T(fvc::grad(U)))
          + ((this->rho_*lambda_)*fvc::div(U))
           *dimensioned<symmTensor>("I", dimless, symmTensor::I)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::dissipationCoeff() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            12.0*(1.0 - sqr(e_))*gs0_*sqr(*this)*rho_
           /(sqrt(Foam::constant::mathematical::pi)*d())
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::productionCoeff() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "dissipationCoeff",
            81.0*(*this)*sqr(otherPhase().mu())*magSqr(U_ - otherPhase().U())
           /(
                gs0_*pow3(d())*rho_*sqrt(Foam::constant::mathematical::pi)
              + dimensionedScalar("small", dimMass, small)
            )
        )
    );
}


void Foam::granularPhaseModel::advect
(
    const dimensionedScalar& deltaT,
    const dimensionedVector& g,
    const volVectorField& Ui,
    const volScalarField& pi,
    const bool oldTime
)
{
    volScalarField& alpha = *this;
    if (oldTime)
    {
        alphaRho_ = alphaRho_.oldTime() - deltaT*fvc::div(massFlux_);
        alphaRho_.correctBoundaryConditions();

        alphaRhoU_ =
            alphaRhoU_.oldTime()
          - deltaT
           *(
                fvc::div(momentumFlux_)
              + alpha*fvc::grad(pi)
              - alphaRho_*g
            );
        alphaRhoU_.correctBoundaryConditions();

        alphaRhoE_ =
            alphaRhoE_.oldTime()
          - deltaT*fvc::div(energyFlux_);
        alphaRhoE_.correctBoundaryConditions();

        alphaRhoPTE_ =
            alphaRhoPTE_.oldTime()
          - deltaT*(fvc::div(PTEFlux_)  + Ps_*fvc::div(phiPtr_()));
        alphaRhoPTE_.correctBoundaryConditions();
    }
    else
    {
        alphaRho_ -= deltaT*fvc::div(massFlux_);
        alphaRho_.correctBoundaryConditions();

        alphaRhoU_ -=
            deltaT
           *(
                fvc::div(momentumFlux_)
              + alpha*fvc::grad(pi)
              - alphaRho_*g
            );
        alphaRhoU_.correctBoundaryConditions();

        alphaRhoE_ -= deltaT*(fvc::div(energyFlux_));
        alphaRhoE_.correctBoundaryConditions();

        alphaRhoPTE_ -= deltaT*(fvc::div(PTEFlux_) + Ps_*fvc::div(phiPtr_()));
        alphaRhoPTE_.correctBoundaryConditions();
    }
}


void Foam::granularPhaseModel::solveSources
(
    const dimensionedScalar& deltaT,
    const volScalarField& Mdot,
    const volVectorField& Ftot,
    const volScalarField& qConv,
    const volScalarField& gamma,
    const volVectorField& Ui,
    const volScalarField& pi
)
{
    volScalarField dissipation
    (
        (12.0*(1.0 - sqr(e_))*gs0_*sqr(*this)*rho_*pow(Theta_, 1.5))
       /(sqrt(Foam::constant::mathematical::pi)*d())
    );

    alphaRho_ += deltaT*Mdot;
    alphaRhoU_ += deltaT*Ftot;
    alphaRhoE_ += deltaT*(dissipation + qConv);
    alphaRhoPTE_ += deltaT*(-dissipation + gamma);
}


void Foam::granularPhaseModel::updateFluxes()
{
    // calculate fluxes with
    volScalarField Ptot
    (
        IOobject::groupName("Ptot", name()),
        Ps_ + Pfric_
    );
    this->fluxFunction_->updateFluxes
    (
        alphaf_,
        massFlux_,
        momentumFlux_,
        energyFlux_,
        PTEFlux_,
        *this,
        rho_,
        U_,
        E_,
        Theta_,
        Ptot,
        c(),
        fluid_.p()
    );
    gradAlpha_ = fvc::surfaceIntegrate(fluid_.mesh().Sf()*alphaf_);
}


void Foam::granularPhaseModel::decode()
{
    volScalarField& alpha = *this;
    alpha = alphaRho_/rho_;
    alpha.max(residualAlpha_);
    alphaRho_ = alpha*rho_;

    U_ = alphaRhoU_/alphaRho_;
    phiPtr_() = fvc::flux(U_);

    E_ = alphaRhoE_/this->alphaRho_;
    he_ = E_;
    thermoPtr_->correct();

    // Update kinetic theory quantities
    Theta_ = 2.0/3.0*alphaRhoPTE_/alphaRho_;
    Theta_.max(0);
    Theta_.min(100);

    volScalarField ThetaSqrt(sqrt(Theta_));
    scalar sqrtPi = sqrt(Foam::constant::mathematical::pi);
    volScalarField da(dPtr_->d());

    gs0_ = radialModel_->g0(*this, alphaMinFriction_, alphaMax_);
    kappa_ = conductivityModel_->kappa(*this, Theta_, gs0_, rho_, da, e_);
    nut_ = viscosityModel_->nu(*this, Theta_, gs0_, rho_, da, e_);
    lambda_ = (4.0/3.0)*sqr(alpha)*da*gs0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

    // Frictional pressure
    Pfric_ =
        frictionalStressModel_->frictionalPressure
        (
            *this,
            alphaMinFriction_,
            alphaMax_
        );

    volSymmTensorField D(symm(fvc::grad(U_)));
    nuFric_ = frictionalStressModel_->nu
    (
        *this,
        alphaMinFriction_,
        alphaMax_,
        Pfric_/rho_,
        D
    );

    // Limit viscosity and add frictional viscosity
    nut_.min(maxNut_);
    nuFric_ = Foam::min(nuFric_, maxNut_ - nut_);
    nut_ += nuFric_;

    Ps_ =
        granularPressureModel_->granularPressureCoeff
        (
            *this,
            gs0_,
            rho_,
            e_
        )*Theta_;
}


void Foam::granularPhaseModel::encode()
{
    (*this).max(residualAlpha_);
    (*this).correctBoundaryConditions();

    alphaRho_ = (*this)*rho_;
    alphaRho_.correctBoundaryConditions();

    alphaRhoU_ = alphaRho_*U_;
    alphaRhoU_.correctBoundaryConditions();

    E_ = he_;
    alphaRhoE_ = alphaRho_*E_;
    alphaRhoE_.correctBoundaryConditions();

    alphaRhoPTE_ = 1.5*alphaRho_*Theta_;
    alphaRhoPTE_.correctBoundaryConditions();
}


void Foam::granularPhaseModel::correctThermo()
{
    // Local references
    volScalarField alpha(*this);
    alpha.max(residualAlpha_);

    const volScalarField& rho = this->rho();
    const volVectorField& U = U_;
    const volVectorField& Uc = otherPhase().U();

    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    dimensionedScalar ThetaSmall("ThetaSmall", Theta_.dimensions(), 1.0e-6);
    dimensionedScalar ThetaSmallSqrt(sqrt(ThetaSmall));

    tmp<volScalarField> tda(d());
    const volScalarField& da = tda();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));

    // Calculating the radial distribution function
    gs0_ = radialModel_->g0(alpha, alphaMinFriction_, alphaMax_);

    // Particle viscosity (Table 3.2, p.47)
    nut_ = viscosityModel_->nu(alpha, Theta_, gs0_, rho, da, e_);

    volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));

    // Bulk viscosity  p. 45 (Lun et al. 1984).
    lambda_ = (4.0/3.0)*sqr(alpha)*da*gs0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

    // Stress tensor, Definitions, Table 3.1, p. 43
    volSymmTensorField tau
    (
        rho*(2.0*nut_*D + (lambda_ - (2.0/3.0)*nut_)*tr(D)*I)
    );

    // Dissipation (Eq. 3.24, p.50)
    volScalarField gammaCoeff
    (
        "gammaCoeff",
        12.0*(1.0 - sqr(e_))
        *Foam::max(sqr(alpha), residualAlpha_)
        *rho*gs0_*(1.0/da)*ThetaSqrt/sqrtPi
    );

    // Drag
    volScalarField beta
    (
        fluid_.lookupSubModel<dragModel>(*this, otherPhase()).K()
    );

    // Eq. 3.25, p. 50 Js = J1 - J2
    volScalarField J1("J1", 3.0*beta);
    volScalarField J2
    (
        "J2",
        0.25*sqr(beta)*da*magSqr(U - Uc)
        /(
            Foam::max(alpha, residualAlpha_)*rho
            *sqrtPi*(ThetaSqrt + ThetaSmallSqrt)
        )
    );

    // particle pressure - coefficient in front of Theta (Eq. 3.22, p. 45)
    volScalarField PsCoeff
    (
        granularPressureModel_->granularPressureCoeff
        (
            alpha,
            gs0_,
            rho,
            e_
        )
    );

    // 'thermal' conductivity (Table 3.3, p. 49)
    kappa_ = conductivityModel_->kappa(alpha, Theta_, gs0_, rho, da, e_);


    // Construct the granular temperature equation (Eq. 3.20, p. 44)
    // NB. note that there are two typos in Eq. 3.20:
    //     Ps should be without grad
    //     the laplacian has the wrong sign
//     fvScalarMatrix ThetaEqn
//     (
//         1.5*
//         (
//             fvm::ddt(alpha, rho, Theta_)
//             + fvm::div(alphaRhoPhi, Theta_)
//             - fvc::Sp(fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi), Theta_)
//         )
//         - fvm::laplacian(kappa_, Theta_, "laplacian(kappa,Theta)")
//         ==
//         - fvm::SuSp((PsCoeff*I) && gradU, Theta_)
//         + (tau && gradU)
//         + fvm::Sp(-gammaCoeff, Theta_)
//         + fvm::Sp(-J1, Theta_)
//         + fvm::Sp(J2/(Theta_ + ThetaSmall), Theta_)
//         + fvOptions(alpha, rho, Theta_)
//     );


    Theta_.max(0);
    Theta_.min(100);

    // particle viscosity (Table 3.2, p.47)
    nut_ = viscosityModel_->nu(alpha, Theta_, gs0_, rho, da, e_);

    ThetaSqrt = sqrt(Theta_);

    // Bulk viscosity  p. 45 (Lun et al. 1984).
    lambda_ = (4.0/3.0)*sqr(alpha)*da*gs0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

    // Frictional pressure
    volScalarField pf
    (
        frictionalStressModel_->frictionalPressure
        (
            *this,
            alphaMinFriction_,
            alphaMax_
        )
    );

    nuFric_ = frictionalStressModel_->nu
    (
        *this,
        alphaMinFriction_,
        alphaMax_,
        pf/rho,
        D
    );

    // Limit viscosity and add frictional viscosity
    nut_.min(maxNut_);
    nuFric_ = Foam::min(nuFric_, maxNut_ - nut_);
    nut_ += nuFric_;

    Ps_ = PsCoeff*Theta_;
    Pfric_ =  pf;

    if (debug)
    {
        Info<< typeName << ':' << nl
            << "    max(Theta) = " << Foam::max(Theta_).value() << nl
            << "    max(nut) = " << Foam::max(nut_).value() << endl;
    }
}


bool Foam::granularPhaseModel::read(const dictionary& phaseProperties)
{
    phaseDict_.lookup("equilibrium") >> equilibrium_;
    e_.readIfPresent(phaseDict_);
    alphaMax_.readIfPresent(phaseDict_);
    alphaMinFriction_.readIfPresent(phaseDict_);

    viscosityModel_->read();
    conductivityModel_->read();
    radialModel_->read();
    granularPressureModel_->read();
    frictionalStressModel_->read();

    return true;
}


void Foam::granularPhaseModel::store()
{
    (*this).storeOldTime();
    alphaRho_.storeOldTime();
    alphaRhoU_.storeOldTime();
    alphaRhoE_.storeOldTime();
    alphaRhoPTE_.storeOldTime();
    Theta_.storeOldTime();
}

// ************************************************************************* //
