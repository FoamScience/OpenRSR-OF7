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
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Only used internally
#define makeTypePhasesTypeName(typePhase, dataType)                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(typePhase<dataType>, 0)

// Sometimes used externally
#define makePhasesTypeName(typePhase)                                          \
                                                                               \
    makeTypePhasesTypeName(typePhase, UnifromMu);                              \
    makeTypePhasesTypeName(typePhase, ChangingMu);

#define makeIncompressiblePhaseType(typePhase, dataType, ns)                   \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        ns::typePhase<dataType>, 0                                             \
    );                                                                         \
    namespace ns {                                                             \
    addTemplatedToRunTimeSelectionTable                                        \
    (                                                                          \
        IncompressiblePhase, typePhase, dataType, dictionary                   \
    );                                                                         \
    }

#define makeCompressiblePhaseType(typePhase, dataType, ns)                     \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        ns::typePhase<dataType>, 0                                             \
    );                                                                         \
    namespace ns {                                                             \
    addTemplatedToRunTimeSelectionTable                                        \
    (                                                                          \
        CompressiblePhase, typePhase, dataType, dictionary                     \
    );                                                                         \
    }

#endif

// ************************************************************************* //
