/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::ReactingSolidPhaseModel

Description
    Class which represents phases with volumetric reactions. Returns the
    reaction rate and heat.

SourceFiles
    ReactingSolidPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingSolidPhaseModel_H
#define ReactingSolidPhaseModel_H

#include "phaseModel.H"
#include "basicSolidChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ReactingSolidPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel, class SolidReactionType, class GasReactionType>
class ReactingSolidPhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- Reaction model
        autoPtr<basicSolidChemistryModel> chemistryPtr_;

        //- Gas phase name
        const word gasPhaseName_;

        //- Reference to gas phase thermo
        const GasReactionType& gasThermo_;


public:

    // Constructors

        ReactingSolidPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );


    //- Destructor
    virtual ~ReactingSolidPhaseModel();


    // Member Functions

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Return the species fraction equation
        virtual tmp<fvScalarMatrix> R(volScalarField& Yi) const;

        //- Return heat release rate
        virtual tmp<volScalarField> Qdot(const bool = true) const;

        //- Return the species fraction equation
        virtual tmp<volScalarField> dmdt(const bool local = true) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingSolidPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
