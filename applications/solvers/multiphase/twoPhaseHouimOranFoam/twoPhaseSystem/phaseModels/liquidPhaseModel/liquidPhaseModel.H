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

Description
    Phase model for a slave pressure, based on a 6 equation model

Class
    Foam::liquidPhaseModel

SourceFiles
    liquidPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef liquidPhaseModel_H
#define liquidPhaseModel_H

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class liquidPhaseModel
:
    public phaseModel
{
    // Private member data

        //- Constant reference to thermodynamic pressure
        const volScalarField& p_;


public:

    // Constructors
        liquidPhaseModel
        (
            const twoPhaseSystem& fluid,
            const dictionary& dict,
            const word& phaseName
        );

    //- Destructor
    virtual ~liquidPhaseModel();


    // Member Functions

        //- Return that the phase uses a slave pressure
        virtual bool slavePressure() const
        {
            return true;
        }

        //- Is the phase granular
        virtual bool granular() const
        {
            return false;
        }

        virtual const volScalarField& p() const
        {
            return p_;
        }

        //- Non-const access to the phase pressure
        virtual volScalarField& p()
        {
            NotImplemented;
            return *this;
        }

        //- Advect with conservative fluxes
        virtual void advect
        (
            const dimensionedScalar& deltaT,
            const dimensionedVector& g,
            const volVectorField& Ui,
            const volScalarField& pi,
            const bool oldTime
        );

        //- Update conservative variables with sources
        virtual void solveSources
        (
            const dimensionedScalar& deltaT,
            const volScalarField& Mdot,
            const volVectorField& Ftot,
            const volScalarField& qConv,
            const volScalarField& gamma,
            const volVectorField& Ui,
            const volScalarField& pi
        );

        //- Decode conservative variables
        virtual void updateFluxes();

        //- Correct primatives based on internal energy/enthalpy
        virtual void correctThermo();

        //- Encode primative variables
        virtual void encode();

        //- Decode conservative variables
        virtual void decode();

        //- Set old time values to current
        virtual void store();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
