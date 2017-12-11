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

Class
    Foam::polydisperseKineticTheoryModel

Description
    Kinetic theory particle phase RAS model

    Reference:
    \verbatim
        "Derivation, implementation, and validation of computer simulation
         models for gas-solid fluidized beds",
        van Wachem, B.G.M.,
        Ph.D. Thesis, Delft University of Technology, Amsterdam, 2000.
    \endverbatim

    There are no default model coefficients.

SourceFiles
    polydisperseKineticTheoryModel.C

\*---------------------------------------------------------------------------*/

#ifndef polydisperseKineticTheoryModel_H
#define polydisperseKineticTheoryModel_H

#include "phaseSystem.H"
#include "phaseModel.H"
#include "phasePair.H"
#include "HashPtrTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace RASModels
{
    class kineticTheoryModel;
}

namespace kineticTheoryModels
{
    class radialModel;
    class frictionalStressModel;
//     class granularPressureModel;
}

/*---------------------------------------------------------------------------*\
                 Class polydisperseKineticTheoryModel Declaration
\*---------------------------------------------------------------------------*/

class polydisperseKineticTheoryModel
:
    public regIOobject
{
    // Private data

        //- Const refrence to multiphaseSystem
        const phaseSystem& fluid_;

        //- Kinetic theory dictionary
        const dictionary& dict_;

        //- Name of total partical phase
        const word name_;

        //- List of granular phase names
        wordList phases_;

        //- Table of granular pairs
        HashTable
        <
            autoPtr<phasePair>,
            phasePairKey,
            phasePairKey::hash
        > pairs_;

        //- Run-time selected radial distribution model
        autoPtr<kineticTheoryModels::radialModel> radialModel_;

        //- Run-time selected frictional stress model
        autoPtr<kineticTheoryModels::frictionalStressModel>
            frictionalStressModel_;


        // Kinetic Theory Model coefficients

            //- Table for coefficient of restitutions based on phase pairs
            phasePair::scalarTable eTable_;

            //- Table for coefficient of fricitions based on phase pairs
            phasePair::scalarTable CfTable_;

            //- Maximum packing phase-fraction
            volScalarField alphaMax_;

            //- Residual phase fraction
            dimensionedScalar residualAlpha_;


        // Total granualr fields

            //- Total granular phase-fraction
            volScalarField alphap_;

            //- Phase averaged velocity
            volVectorField Up_;


        // Stored tables of coefficients

            //- Granular pressure coeffients for granular pairs
            HashPtrTable
            <
                volScalarField,
                phasePairKey,
                phasePairKey::hash
            > PsCoeffs_;

            //- Thermal conductivity between granular pairs
            HashPtrTable
            <
                volScalarField,
                phasePairKey,
                phasePairKey::hash
            > kappas_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        polydisperseKineticTheoryModel
        (
            const polydisperseKineticTheoryModel& fluid
        );

        //- Disallow default bitwise assignment
        void operator=(const polydisperseKineticTheoryModel&);

        //- Compute the packing limit for phase1 and phase2
        scalar calcAlphaMax(const label, const label, const scalarList&) const;


public:

//     //- Runtime type information
//     TypeName("polydisperseKineticTheory");


    // Constructors

        //- Construct from components
        polydisperseKineticTheoryModel(const phaseSystem&);


    //- Destructor
    ~polydisperseKineticTheoryModel();


    // Member Functions

        //- Re-read model coefficients if they have changed
        bool read();

        //- Return fluid
        const phaseSystem& fluid() const
        {
            return fluid_;
        }

        //- Const refrence to list of granular phases
        const wordList& phases() const;

        //- Add phase to list of kinetic theory phases.
        void addPhase(const RASModels::kineticTheoryModel&);

        //- Return true if phase is granular
        bool found(const word&) const;

        //- Return total granular phase volume fraction
        const volScalarField& alphap() const
        {
            return alphap_;
        }

        //- Return phase averaged velocity
        const volVectorField& Up() const
        {
            return Up_;
        }

        //- Return maximum volume fraction
        const volScalarField& alphaMax() const
        {
            return alphaMax_;
        }

        //- Coefficients of restitutions between granular phases
        const phasePair::scalarTable es() const
        {
            return eTable_;
        }

        //- coefficients of friction between granular phases
        const phasePair::scalarTable Cfs() const
        {
            return CfTable_;
        }

        //- Return radial distribution function
        tmp<volScalarField> gs0(const phaseModel&, const phaseModel&) const;

        //- Return radial distribution function
        tmp<volScalarField> gs0Prime
        (
            const phaseModel&,
            const phaseModel&
        ) const;

        //- Return the sum of all granular pressure for a given phase
        tmp<volScalarField> PsCoeff(const phaseModel&) const;

        //- Return the sum of the derivative of granular pressure
        //  for a given phase
        tmp<volScalarField> PsCoeffPrime(const phaseModel&) const;

        //- Return the fricitional-pressure
        tmp<volScalarField> frictionalPressure(const phaseModel&) const;

        //- Return the frictional-pressure'
        // (derivative of frictional-pressure w.r.t. phase-fraction)
        tmp<volScalarField> frictionalPressurePrime(const phaseModel&) const;

        //- Return the frictional-viscosity'
        tmp<volScalarField> nuFrictional(const phaseModel&) const;

        //- Solve the kinetic theory equations and correct the viscosity
        void correct();

        //- Correct volume fraction total
        void correctAlphap();

        virtual bool writeData(Ostream&) const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
