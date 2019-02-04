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

Class
    Foam::dragOde

Description
    Basic ode solver for coupling two phases with drag

SourceFiles
    dragOdeI.H
    dragOde.C

\*---------------------------------------------------------------------------*/

#ifndef dragOde_H
#define dragOde_H

#include "ODESystem.H"
#include "ODESolver.H"
#include "volFields.H"
#include "simpleMatrix.H"
#include "phaseSystem.H"
#include "labelVector.H"
#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class dragOde Declaration
\*---------------------------------------------------------------------------*/

class dragOde
:
    public ODESystem
{
    typedef PtrListDictionary<phaseModel> phaseModelList;
    typedef HashTable
        <
            autoPtr<phasePair>,
            phasePairKey,
            phasePairKey::hash
        > phasePairTable;
    typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<dragModel>>,
            phasePairKey,
            phasePairKey::hash
        > dragModelTable;

    // Private Member Functions

        //- Update velocities in phase models
        void setUs();

        //- Update q vector
        void setq();

        //- Disallow copy constructor
        dragOde(const dragOde&);

        //- Disallow default bitwise assignment
        void operator=(const dragOde&);


protected:


    // Protected data

        //- Ode dictionary
        const dictionary& dict_;

        //- Solve drag with ode
        Switch solveDrag_;

        //- Reference to phase system
        const phaseSystem& fluid_;

        //- Reference to phase models
        phaseModelList& phaseModels_;

        //- PrtList of references to dragModels
        PtrList<BlendedInterfacialModel<dragModel>> dragModels_;

        //- Reference to phase pairs
        PtrList<phasePair> phasePairs_;

        //- vector of solution dimensions
        labelVector solutionD_;

        //- number of solution dimensions
        label nDims_;

        //- Number of equations solved
        label nEqns_;

        //- Temporary field (u1 v1 w1, u2 v2 w2 ...)
        scalarField q_;

        //- Temporary rate-of-change of variables field
        scalarField dqdt_;

        //- Stored time step
        scalarField deltaT_;

        //- Current cell
        label celli_ = 0;

        //- Ode solver
        autoPtr<ODESolver> odeSolver_;


public:


    // Constructors

        //- Construct from thermo
        dragOde(phaseSystem& fluid, dragModelTable& dragModels);


    //- Destructor
    virtual ~dragOde();


    // Member Functions

        bool solveDrag() const
        {
            return solveDrag_;
        }

        //- Solve the reaction system for the given time step
        //  and return the characteristic time
        virtual scalar solve(const scalar& deltaT);

    // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            virtual label nEqns() const
            {
                return nEqns_;
            }

            virtual void derivatives
            (
                const scalar t,
                const scalarField& q,
                scalarField& dqdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& q,
                scalarField& dqdt,
                scalarSquareMatrix& J
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
