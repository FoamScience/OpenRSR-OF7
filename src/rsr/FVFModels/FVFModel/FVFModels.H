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
    FVFModels.C

\*---------------------------------------------------------------------------*/

#ifndef FVFModels_H
#define FVFModels_H

#include "FVFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Only used internally
#define makeTypeFVFModelsTypeName(typeFVFModel, dataType)                      \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(typeFVFModel<dataType>, 0)

// Sometimes used externally
#define makeFVFModelsTypeName(typeFVFModel)                                    \
                                                                               \
    makeTypeFVFModelsTypeName(typeFVFModel, Incompressible);                   \
    makeTypeFVFModelsTypeName(typeFVFModel, Compressible);                     \

// Define type info for single dataType template instantiation
#define makeFVFModelType(typeFVFModel, dataType, ns)                           \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        ns::typeFVFModel<dataType>, 0                                          \
    );                                                                         \
    namespace ns {                                                             \
    addTemplatedToRunTimeSelectionTable                                        \
    (                                                                          \
        FVFModel, typeFVFModel, dataType, dictionary                           \
    );                                                                         \
    }

// Define type info for Compressible, Incompressible etc. instantiations
#define makeFVFModels(typeFVFModel, ns)                                        \
                                                                               \
    makeFVFModelType(typeFVFModel, Compressible, ns);                          \
    makeFVFModelType(typeFVFModel, Incompressible, ns);                             

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
