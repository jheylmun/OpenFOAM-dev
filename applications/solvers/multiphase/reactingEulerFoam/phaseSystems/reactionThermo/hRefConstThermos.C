/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 OpenFOAM Foundation
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

#include "makeReactionThermo.H"
#include "makeSolidThermo.H"
#include "makeReactingSolidThermo.H"
#include "makeThermo.H"
#include "makeChemistryModel.H"
#include "makeChemistrySolverTypes.H"
#include "makeChemistryTabulationMethods.H"
#include "makeChemistryReductionMethods.H"

#include "rhoReactionThermo.H"
#include "solidReactionThermo.H"
#include "heRhoThermo.H"
#include "heSolidThermo.H"
#include "thermo.H"
#include "solidChemistryModel.H"

#include "specie.H"
#include "perfectGas.H"
#include "perfectFluid.H"
#include "rhoConst.H"
#include "hConstThermo.H"
#include "hPowerThermo.H"
#include "constIsoSolidTransport.H"
#include "constAnIsoSolidTransport.H"
#include "exponentialSolidTransport.H"

#include "sensibleEnthalpy.H"

#include "hRefConstThermo.H"

#include "constTransport.H"

#include "pureMixture.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"

#include "thermoPhysicsTypes.H"
#include "solidThermoPhysicsTypes.H"

#include "StandardChemistryModel.H"
#include "TDACChemistryModel.H"
#include "thermoPhysicsTypes.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Thermo type typedefs:

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            perfectGas<specie>
        >,
        sensibleEnthalpy
    >
> constRefGasHThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            perfectFluid<specie>
        >,
        sensibleEnthalpy
    >
> constRefFluidHThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            perfectGas<specie>
        >,
        sensibleInternalEnergy
    >
> constRefGasEThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            perfectFluid<specie>
        >,
        sensibleInternalEnergy
    >
> constRefFluidEThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            rhoConst<specie>
        >,
        sensibleInternalEnergy
    >
> constRefRhoConstEThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            rhoConst<specie>
        >,
        sensibleEnthalpy
    >
> constRefRhoConstHThermoPhysics;

// pureMixture, sensibleEnthalpy:

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    perfectFluid,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    rhoConst,
    specie
);


// pureMixture, sensibleInternalEnergy:

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hRefConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hRefConstThermo,
    perfectFluid,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hRefConstThermo,
    rhoConst,
    specie
);


// multiComponentMixture, sensibleInternalEnergy:

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefGasEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefFluidEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefRhoConstEThermoPhysics
);


// multiComponentMixture, sensibleEnthalpy:

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefRhoConstHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefFluidHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefGasHThermoPhysics
);


// Solid thermos
makeSolidThermoPhysicsType
(
    solidThermo,
    heSolidThermo,
    pureMixture,
    hConstSolidThermoPhysics
);

makeSolidThermoPhysicsType
(
    solidThermo,
    heSolidThermo,
    pureMixture,
    hPowerSolidThermoPhysics
);

makeSolidThermoPhysicsType
(
    solidThermo,
    heSolidThermo,
    pureMixture,
    hTransportThermoPoly8SolidThermoPhysics
);

makeSolidThermoPhysicsType
(
    solidThermo,
    heSolidThermo,
    pureMixture,
    hExpKappaConstSolidThermoPhysics
);

// Solid reaction thermo
makeReactingSolidThermo
(
    solidReactionThermo,
    heSolidThermo,
    reactingMixture,
    constIsoSolidTransport,
    sensibleEnthalpy,
    hConstThermo,
    rhoConst,
    specie
);

makeReactingSolidThermo
(
    solidReactionThermo,
    heSolidThermo,
    reactingMixture,
    constIsoSolidTransport,
    sensibleEnthalpy,
    hPowerThermo,
    rhoConst,
    specie
);

makeReactingSolidThermo
(
    solidThermo,
    heSolidThermo,
    multiComponentMixture,
    constIsoSolidTransport,
    sensibleEnthalpy,
    hConstThermo,
    rhoConst,
    specie
);

// Solid chemisty models
makeChemistryModel(solidReactionThermo);

makeChemistryModelType
(
    StandardChemistryModel,
    solidReactionThermo,
    hConstSolidThermoPhysics
);

makeChemistryModelType
(
    StandardChemistryModel,
    solidReactionThermo,
    hPowerSolidThermoPhysics
);

makeChemistryModelType
(
    StandardChemistryModel,
    solidReactionThermo,
    hTransportThermoPoly8SolidThermoPhysics
);

makeChemistryModelType
(
    StandardChemistryModel,
    solidReactionThermo,
    hExpKappaConstSolidThermoPhysics
);

// Solid chemistry solvers
makeStandardChemistrySolverTypes
(
    solidReactionThermo,
    hConstSolidThermoPhysics
);
makeStandardChemistrySolverTypes
(
    solidReactionThermo,
    hPowerSolidThermoPhysics
);
makeStandardChemistrySolverTypes
(
    solidReactionThermo,
    hTransportThermoPoly8SolidThermoPhysics
);
makeStandardChemistrySolverTypes
(
    solidReactionThermo,
    hExpKappaConstSolidThermoPhysics
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
