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
    basicInterpolationTables.C

\*---------------------------------------------------------------------------*/

#ifndef basicInterpolationTables_H
#define basicInterpolationTables_H

#include "basicInterpolationTable.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Only used internally
#define makeTypeInterpolationTablesTypeName(typeInterpolationTable, dataType)  \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(typeInterpolationTable<dataType>, 0)

// Sometimes used externally
#define makeInterpolationTablesTypeName(typeInterpolationTable)                \
                                                                               \
    makeTypeInterpolationTablesTypeName(typeInterpolationTable, scalar);       \
    makeTypeInterpolationTablesTypeName(typeInterpolationTable, vector);       \
    makeTypeInterpolationTablesTypeName                                        \
    (                                                                          \
        typeInterpolationTable, sphericalTensor                                \
    );                                                                         \
    makeTypeInterpolationTablesTypeName(typeInterpolationTable, symmTensor);   \
    makeTypeInterpolationTablesTypeName(typeInterpolationTable, tensor)

// Define type info for single dataType template instantiation (eg, vector)
#define makeInterpolationTableType(typeInterpolationTable, dataType, ns)       \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        ns::typeInterpolationTable<dataType>, 0                                \
    );                                                                         \
    namespace ns {                                                             \
    addTemplatedToRunTimeSelectionTable                                        \
    (                                                                          \
        basicInterpolationTable, typeInterpolationTable, dataType, dictionary  \
    );                                                                         \
    }

// Define type info for scalar, vector etc. instantiations
#define makeInterpolationTables(typeInterpolationTable, ns)                    \
                                                                               \
    makeInterpolationTableType(typeInterpolationTable, scalar, ns);            \
    makeInterpolationTableType(typeInterpolationTable, vector, ns);            \
    makeInterpolationTableType(typeInterpolationTable, sphericalTensor, ns);   \
    makeInterpolationTableType(typeInterpolationTable, symmTensor, ns);        \
    makeInterpolationTableType(typeInterpolationTable, tensor, ns)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
