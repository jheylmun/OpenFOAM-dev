/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
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
    Foam::stiffenedEoS

Description
    Adiabatic perfect fluid equation of state.

SourceFiles
    stiffenedEoSI.H
    stiffenedEoS.C

\*---------------------------------------------------------------------------*/

#ifndef stiffenedEoS_H
#define stiffenedEoS_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class stiffenedEoS;

template<class Specie>
inline stiffenedEoS<Specie> operator+
(
    const stiffenedEoS<Specie>&,
    const stiffenedEoS<Specie>&
);

template<class Specie>
inline stiffenedEoS<Specie> operator*
(
    const scalar,
    const stiffenedEoS<Specie>&
);

template<class Specie>
inline stiffenedEoS<Specie> operator==
(
    const stiffenedEoS<Specie>&,
    const stiffenedEoS<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const stiffenedEoS<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class stiffenedEoS Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class stiffenedEoS
:
    public Specie
{
    // Private data

        //- Reference pressure
        scalar p0_;

        //- Reference density
        scalar rho0_;

        //- The isentropic exponent
        scalar gamma_;

        //- Pressure offset for a stiffened gas
        scalar B_;


public:

    // Constructors

        //- Construct from components
        inline stiffenedEoS
        (
            const Specie& sp,
            const scalar p0,
            const scalar rho0,
            const scalar gamma,
            const scalar B
        );

        //- Construct from dictionary
        stiffenedEoS(const dictionary& dict);

        //- Construct as named copy
        inline stiffenedEoS
        (
            const word& name,
            const stiffenedEoS&
        );

        //- Construct and return a clone
        inline autoPtr<stiffenedEoS> clone() const;

        // Selector from dictionary
        inline static autoPtr<stiffenedEoS> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "stiffenedEoS<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- return pressure [Pa]
            inline scalar p(scalar T, scalar rho) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;

            //- Return the derivative of pressure with respesct to temperature
            inline scalar dpdT
            (
                scalar T,
                scalar rho
            ) const;

            //- Return the derivative of pressure with respesct to molar volume
            inline scalar dpdv
            (
                scalar T,
                scalar rho
            ) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const stiffenedEoS&);
        inline void operator*=(const scalar);


    // Friend operators

        friend stiffenedEoS operator+ <Specie>
        (
            const stiffenedEoS&,
            const stiffenedEoS&
        );

        friend stiffenedEoS operator* <Specie>
        (
            const scalar s,
            const stiffenedEoS&
        );

        friend stiffenedEoS operator== <Specie>
        (
            const stiffenedEoS&,
            const stiffenedEoS&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const stiffenedEoS&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "stiffenedEoSI.H"

#ifdef NoRepository
    #include "stiffenedEoS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
