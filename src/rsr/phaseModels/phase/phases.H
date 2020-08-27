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

SourceFiles
    phases.C

\*---------------------------------------------------------------------------*/

#ifndef phases_H
#define phases_H

#include "UniformityTypes.H"
#include "debug.H"
#include "phase.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

    // A type for constant-viscosity phases
    template<class CompressibilityType>
    using ConstantMuPhase = phase<CompressibilityType, UniformMu>;

    // A type for variant-viscosity phases
    template<class CompressibilityType>
    using ChangingMuPhase = phase<CompressibilityType, ChangingMu>;
}

// Only used internally
#define makeTypePhasesTypeName(typePhase, dataTypes...)                        \
                                                                               \
    template<>                                                                 \
    const ::Foam::word typePhase<dataTypes>::typeName                          \
        ( typePhase<dataTypes>::typeName_() );                                 \
    template<>                                                                 \
    int typePhase<dataTypes>::debug                                            \
    (                                                                          \
        ::Foam::debug::debugSwitch                                             \
        (                                                                      \
            typePhase<dataTypes>::typeName_(), 0                               \
        )                                                                      \
    );                                                                         \
    template<>                                                                 \
    const Foam::RegisterDebugSwitch<typePhase<dataTypes>>                      \
        Foam::RegisterDebugSwitch<typePhase<dataTypes>>::registerDebugSwitch   \
        ( typePhase<dataTypes>::typeName_() )

// Sometimes used externally
#define makePhasesTypeName(typePhase)                                          \
                                                                               \
    makeTypePhasesTypeName(typePhase, Compressible, UniformMu);                \
    makeTypePhasesTypeName(typePhase, Incompressible, UniformMu);              \
    makeTypePhasesTypeName(typePhase, Compressible, ChangingMu)

// Define type info for a template instantiation
#define makeConstantMuPhaseType(typePhase, dataType, ns)                       \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        ns::typePhase<dataType>, 0                                             \
    );                                                                         \
    namespace ns {                                                             \
    addTemplatedToRunTimeSelectionTable                                        \
    (                                                                          \
        ConstantMuPhase, typePhase, dataType, dictionary                       \
    );                                                                         \
    }

// Define type info for a template instantiation
#define makeChangingMuPhaseType(typePhase, dataType, ns)                       \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        ns::typePhase<dataType>, 0                                             \
    );                                                                         \
    namespace ns {                                                             \
    addTemplatedToRunTimeSelectionTable                                        \
    (                                                                          \
        ChangingMuPhase, typePhase, dataType, dictionary                       \
    );                                                                         \
    }

// Define type info for scalar, vector etc. instantiations
#define makePhases(typePhase, ns)                    \
                                                                               \
    makeConstantMuPhaseType(typePhase, Compressible, ns);                      \
    makeConstantMuPhaseType(typePhase, Incompressible, ns);                      

#endif

// ************************************************************************* //
