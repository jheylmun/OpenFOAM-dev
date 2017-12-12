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

Class
    Foam::fluxIntegrators::Euler

Description
    1st order Runge-Kutta ode time integrater

SourceFiles
    RK45.C

\*---------------------------------------------------------------------------*/

#ifndef Euler_H
#define Euler_H

#include "fluxIntegrator.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxIntegrators
{

/*---------------------------------------------------------------------------*\
                        Class Euler Declaration
\*---------------------------------------------------------------------------*/

class Euler
:
    public fluxIntegrator
{
    // Private data

public:

    //- Runtime type information
    TypeName("Euler");

    // Constructor
    Euler(compressibleSystem& fluid, const word& phaseName);

    //- Destructor
    virtual ~Euler();


    // Edit functions

        //- Transport conservative variables
        virtual void integrateFluxes
        (
            volScalarField& rho,
            volVectorField& rhoU,
            volScalarField& rhoE
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluxIntegrators
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
