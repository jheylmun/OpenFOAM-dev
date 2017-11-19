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

#include "heRhoThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::calculate()
{
    const scalarField& hCells = this->he();
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& rhoBf =
        this->rho_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heRhoThermo<BasicPsiThermo, MixtureType>::heRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName)
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heRhoThermo<BasicPsiThermo, MixtureType>::~heRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::correctP()
{
    const scalarField& hCells = this->he();
    const scalarField& rhoCells = this->rho_.primitiveFieldRef();

    scalarField& pCells = this->p_;
    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        TCells[celli] = mixture_.THErho
        (
            hCells[celli],
            rhoCells[celli],
            TCells[celli]
        );

        pCells[celli] = mixture_.p(TCells[celli], rhoCells[celli]);
        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& rhoBf =
        this->rho_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pp[facei] = mixture_.p(pT[facei], prho[facei]);
                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.THErho
                (
                    phe[facei],
                    prho[facei],
                    pT[facei]
                );

                pp[facei] = mixture_.p(pT[facei], prho[facei]);
                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }
    this->p_.correctBoundaryConditions();
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::calcRho()
{
    scalarField& rhoInternal(this->rho_.primitiveFieldRef());

    forAll(this->rho_, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);
        rhoInternal[celli] =
            mixture_.rho
            (
                this->p_[celli],
                this->T_[celli]
            );
    }

    const volScalarField::Boundary& TBf = this->T_.boundaryField();
    const volScalarField::Boundary& pBf = this->p_.boundaryField();
    volScalarField::Boundary& rhoBf = this->rho_.boundaryFieldRef();

    forAll(this->rho_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = TBf[patchi];
        const fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];

        forAll(pT, facei)
        {
            const typename MixtureType::thermoType& mixture_ =
                this->patchFaceMixture(patchi, facei);

            prho[facei] = mixture_.rho(pp[facei], pT[facei]);
        }
    }
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::heRho()
{
    const scalarField& pCells = this->p();
    const scalarField& rhoCells = this->rho_.primitiveFieldRef();

    scalarField& heCells = this->he_;
    scalarField& TCells = this->T_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        TCells[celli] = mixture_.TpRho
        (
            pCells[celli],
            rhoCells[celli],
            TCells[celli]
        );

        heCells[celli] = mixture_.HE(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& rhoBf =
        this->rho_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.TpRho
                (
                    pp[facei],
                    prho[facei],
                    pT[facei]
                );
                phe[facei] = mixture_.HE(pp[facei], pT[facei]);
            }
        }
    }
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heRhoThermo<BasicPsiThermo, MixtureType>::dHEdp() const
{
    tmp<volScalarField> tmpdhedp
    (
        new volScalarField
        (
            IOobject
            (
                "dhedp",
                this->p_.time().timeName(),
                this->p_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->p_.mesh(),
            dimensionedScalar
            (
                "0",
                this->he_.dimensions()/dimPressure,
                0.0
            )
        )
    );
    volScalarField& dhedp(tmpdhedp.ref());

    forAll(this->T_, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);
        dhedp[celli] =
            mixture_.dHEdp
            (
                this->T_[celli],
                this->rho_[celli]
            );
    }

    const volScalarField::Boundary& TBf = this->T_.boundaryField();
    const volScalarField::Boundary& rhoBf = this->rho_.boundaryField();
    volScalarField::Boundary& dhedpBf = dhedp.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = TBf[patchi];
        const fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& pdhedp = dhedpBf[patchi];

        forAll(pT, facei)
        {
            const typename MixtureType::thermoType& mixture_ =
                this->patchFaceMixture(patchi, facei);

            pdhedp[facei] = mixture_.dHEdp(pT[facei], prho[facei]);
        }
    }
    return tmpdhedp;
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heRhoThermo<BasicPsiThermo, MixtureType>::dHEdrho() const
{
    tmp<volScalarField> tmpdHEdrho
    (
        new volScalarField
        (
            IOobject
            (
                "dHEdrho",
                this->p_.time().timeName(),
                this->p_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->p_.mesh(),
            dimensionedScalar
            (
                "0",
                this->he_.dimensions()/dimDensity,
                0.0
            )
        )
    );
    volScalarField& dHEdrho(tmpdHEdrho.ref());

    forAll(this->T_, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);
        dHEdrho[celli] =
            mixture_.dHEdrho
            (
                this->T_[celli],
                this->rho_[celli]
            );
    }

    const volScalarField::Boundary& TBf = this->T_.boundaryField();
    const volScalarField::Boundary& rhoBf = this->rho_.boundaryField();
    volScalarField::Boundary& dHEdrhoBf = dHEdrho.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = TBf[patchi];
        const fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& pdHEdrho = dHEdrhoBf[patchi];

        forAll(pT, facei)
        {
            const typename MixtureType::thermoType& mixture_ =
                this->patchFaceMixture(patchi, facei);

            pdHEdrho[facei] = mixture_.dHEdrho(pT[facei], prho[facei]);
        }
    }
    return tmpdHEdrho;
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


// ************************************************************************* //
