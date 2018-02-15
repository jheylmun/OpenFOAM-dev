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
    Foam::fluxIntegrator

Description
    Class used to handle encoding/decoding of primative and conservative
    variables used in compressible flows. Also includes funcitons to calculate
    the conservative fluxes

SourceFiles
    fluxIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef fluxIntegrator_H
#define fluxIntegrator_H

#include "IOdictionary.H"
#include "compressibleSystem.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class fluxIntegrator Declaration
\*---------------------------------------------------------------------------*/

class fluxIntegrator
{
protected:

    // Protected data

        //- Name of compressible phase ('\0' by default)
        compressibleSystem& fluid_;

        //  Conservative fluxes

            //- Mass flux
            surfaceScalarField& massFlux_;

            //- Momentum flux
            surfaceVectorField& momentumFlux_;

            //- Energy flux
            surfaceScalarField& energyFlux_;

public:

    //- Runtime type information
    TypeName("fluxIntegrator");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluxIntegrator,
            dictionary,
            (
                compressibleSystem& fluid,
                const word& phaseName
            ),
            (fluid, phaseName)
        );


    // Constructor
    fluxIntegrator(compressibleSystem& fluid, const word& name);

    // Selector
    static autoPtr<fluxIntegrator> New
    (
        compressibleSystem&,
        const word& phaseName
    );


    //- Destructor
    virtual ~fluxIntegrator();


    // Member Functions

        //- Transport conservative variables
        virtual void integrateFluxes
        (
            volScalarField& rhoOld,
            volVectorField& rhoUOld,
            volScalarField& rhoEOld
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
