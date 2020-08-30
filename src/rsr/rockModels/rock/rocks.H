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

\*---------------------------------------------------------------------------*/

#ifndef rocks_H
#define rocks_H

#include "UniformityTypes.H"
#include "debug.H"
#include "rock.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Only used internally
#define makeTypeRocksTypeName(typeRock, dataType)                              \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(typeRock<dataType>, 0)

// Sometimes used externally
#define makeRocksTypeName(typeRock)                                            \
                                                                               \
    makeTypeRocksTypeName(typeRock, UnifromMu);                                \
    makeTypeRocksTypeName(typeRock, ChangingMu);

#define makeIsoRockType(typeRock, dataType, ns)                                \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        ns::typeRock<dataType>, 0                                              \
    );                                                                         \
    namespace ns {                                                             \
    addTemplatedToRunTimeSelectionTable                                        \
    (                                                                          \
        IsoRock, typeRock, dataType, dictionary                                \
    );                                                                         \
    }

#define makeDiagAnisoRockType(typeRock, dataType, ns)                          \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        ns::typeRock<dataType>, 0                                              \
    );                                                                         \
    namespace ns {                                                             \
    addTemplatedToRunTimeSelectionTable                                        \
    (                                                                          \
        DiagAnisoRock, typeRock, dataType, dictionary                          \
    );                                                                         \
    }

#endif

// ************************************************************************* //
