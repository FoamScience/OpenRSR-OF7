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

#include "scalarList.H"
#include "vectorList.H"
#include "sphericalTensorList.H"
#include "symmTensorList.H"
#include "tensorList.H"

#include "basicInterpolationTables.H"
#include "tableReaders.H"
#include "openFoamTableReader.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define defineTableReaderType(dataType)                                        \
    defineNamedTemplateTypeNameAndDebug(tableReader<dataType >, 0);            \
    defineTemplatedRunTimeSelectionTable(tableReader, dictionary, dataType);

#define makeListTableReaders(typeTableReader)                                  \
                                                                               \
    makeTableReaderType(typeTableReader, scalarList);                          \
    makeTableReaderType(typeTableReader, vectorList);                          \
    makeTableReaderType(typeTableReader, sphericalTensorList);                 \
    makeTableReaderType(typeTableReader, symmTensorList);                      \
    makeTableReaderType(typeTableReader, tensorList)
#define defineInterpolationTableType(dataType)                                 \
    defineNamedTemplateTypeNameAndDebug(basicInterpolationTable<dataType >, 0);\
    defineTemplatedRunTimeSelectionTable                                       \
    (                                                                          \
        basicInterpolationTable,                                               \
        dictionary,                                                            \
        dataType                                                               \
    );

// Define Base Table Readers for list types
defineTableReaderType(scalarList);
defineTableReaderType(vectorList);
defineTableReaderType(sphericalTensorList);
defineTableReaderType(symmTensorList);
defineTableReaderType(tensorList);

// Declare OpenFOAM readers for list tables
makeListTableReaders(openFoamTableReader);

// Define Base Interpolation Tables for primitive types
defineInterpolationTableType(scalar);
defineInterpolationTableType(vector);
defineInterpolationTableType(sphericalTensor);
defineInterpolationTableType(symmTensor);
defineInterpolationTableType(tensor);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
