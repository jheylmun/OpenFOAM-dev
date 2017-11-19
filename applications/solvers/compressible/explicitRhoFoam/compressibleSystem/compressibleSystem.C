/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Jeff Heylmun
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

#include "compressibleSystem.H"
#include "fluxIntegrator.H"
#include "fluxFunction.H"
#include "surfaceInterpolate.H"
#include "fixedValueFvsPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "fvc.H"
#include "fvm.H"


// * * * * * * * * * * * * * * * Protected Functions * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleSystem::compressibleSystem
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    phaseName_(phaseName),
    mesh_(mesh),
    thermoPtr_
    (
        rhoThermo::New(mesh, phaseName)
    ),
    rho_
    (
        IOobject
        (
            IOobject::groupName("rho", phaseName_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermoPtr_().rho(),
        thermoPtr_->T().boundaryField().types()
    ),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", phaseName_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    p_(thermoPtr_->p()),
    E_
    (
        IOobject
        (
            IOobject::groupName("E", phaseName_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermoPtr_->he() + 0.5*magSqr(U_),
        thermoPtr_->he().boundaryField().types()
    ),
    H_
    (
        IOobject
        (
            IOobject::groupName("H", phaseName_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        E_ + p_/rho_
    ),
    rhoU_
    (
        IOobject
        (
            IOobject::groupName("rhoU", phaseName_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*U_,
        U_.boundaryField().types()
    ),
    rhoE_
    (
        IOobject
        (
            IOobject::groupName("rhoE", phaseName_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*E_,
        thermoPtr_->T().boundaryField().types()
    ),
    integrator_
    (
        fluxIntegrator::New
        (
            *this,
            phaseName_
        )
    ),
    fluxFunction_
    (
        fluxFunction::New(mesh_, phaseName_)
    ),
    massFlux_
    (
        IOobject
        (
            IOobject::groupName("massFlux", phaseName_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimVelocity*dimDensity*dimArea, 0)
    ),
    momentumFlux_
    (
        IOobject
        (
            IOobject::groupName("momentumFlux", phaseName_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0", sqr(dimVelocity)*dimDensity*dimArea, Zero)
    ),
    energyFlux_
    (
        IOobject
        (
            IOobject::groupName("energyFlux", phaseName_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", pow3(dimVelocity)*dimDensity*dimArea, 0)
    )
{
    const word phiName = IOobject::groupName("phi", phaseName_);

    IOobject phiHeader
    (
        phiName,
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ
    );

    if (phiHeader.typeHeaderOk<surfaceScalarField>(true))
    {
        Info<< "Reading face flux field " << phiName << endl;

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field " << phiName << endl;

        wordList phiTypes
        (
            U_.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U_.boundaryField(), i)
        {
            if
            (
                isA<fixedValueFvPatchVectorField>(U_.boundaryField()[i])
             || isA<slipFvPatchVectorField>(U_.boundaryField()[i])
             || isA<partialSlipFvPatchVectorField>(U_.boundaryField()[i])
            )
            {
                phiTypes[i] = fixedValueFvsPatchScalarField::typeName;
            }
        }

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("0", dimVelocity*dimArea, 0),//fvc::flux(U_),
                phiTypes
            )
        );
    }
    thermoPtr_->validate("compressibleSystem " + phaseName_, "h", "e");
    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::compressibleSystem::~compressibleSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::compressibleSystem::integrateFluxes()
{
    integrator_->integrateFluxes(rho_, rhoU_, rhoE_);
}


void Foam::compressibleSystem::updateFluxes()
{
    // calculate fluxes with
    fluxFunction_->updateFluxes
    (
        massFlux_,
        momentumFlux_,
        energyFlux_,
        rho_,
        U_,
        H_,
        p_,
        thermoPtr_->speedOfSound()()
    );
}


void Foam::compressibleSystem::decode()
{
    thermoPtr_->rho() = rho_;

    U_ = rhoU_/rho_;
    U_.correctBoundaryConditions();

    phi() = fvc::flux(U_);

    E_ = rhoE_/rho_;

    if (thermoPtr_->he().name()[0] == 'e')
    {
        thermoPtr_->he() = E_ - 0.5*magSqr(U_);
    }
    else if (thermoPtr_->he().name()[0] == 'h')
    {
        NotImplemented
        H_ = E_ + p_/rho_;
        thermoPtr_->he() = H_ - 0.5*magSqr(U_);
    }
    thermoPtr_->correctP();
    H_ = E_ + p_/rho_;
}


void Foam::compressibleSystem::encode()
{
    rho_ = thermoPtr_->rho();
    rho_.correctBoundaryConditions();

    rhoU_ = rho_*U_;
    rhoU_.correctBoundaryConditions();

    rhoE_ = rho_*E_;
    rhoE_.correctBoundaryConditions();
}


void Foam::compressibleSystem::correctThermo()
{
    thermoPtr_->correctP();
    if (thermoPtr_->he().name()[0] == 'e')
    {
        E_ = thermoPtr_->he() + 0.5*magSqr(U_);
        H_ = E_ + p_/rho_;
    }
    else if (thermoPtr_->he().name()[0] == 'h')
    {
        H_ = thermoPtr_->he() + 0.5*magSqr(U_);
        E_ = H_ - p_/rho_;
    }
}


// ************************************************************************* //
