/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "UniformDimensionedField.H"
#include "sourceProperties.H"

namespace Foam
{
    defineTypeNameAndDebug(sourceProperties, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sourceProperties::sourceProperties
(
    const fvMesh& mesh,
    const dictionary& wellDict,
    const cellSet& wellSet,
    const faceSet& faces
)
:
    mesh_(mesh),
    g_
    (
        mesh.lookupObject<UniformDimensionedField<vector>>("g")
    ),
    orientation_
    (
        wordToOrientationHandling
        (
            wellDict.lookup("orientation")
        )
    ),
    operation_
    (
        wordToOpHandling(wellDict.lookup("operationMode"))
    ),
    cells_(wellSet.toc()),
    faces_(faces.toc()),
    V_(wellDict.dictName()+".V", dimVolume, 0.0),
    radius_
    (
        wellDict.dictName()+".radius",
        dimensionedScalar("radius", dimLength, wellDict)), skin_
    (
        readScalar(wellDict.lookup("skin"))
    ),
    J_(),
    injPhase_
    (
     operation_ == operationHandling::injection
     ? word(wellDict.lookup("injectedPhase"))
     : "none"
    )
{
    cellsVolume();
    if (debug and cells_.empty())
    {
        WarningInFunction
            << "No cells in well set " << wellSet.name()
            << ". If a segfault occurs this is the cause." << nl << endl;
    }
}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

typename Foam::sourceProperties::orientationHandling 
Foam::sourceProperties::wordToOrientationHandling
(
    const word& ori
)
{
    if (ori == "vertical")
    {
        return orientationHandling::vertical;
    }
    else if (ori == "horizontalX")
    {
        return orientationHandling::horizontalX;
    }
    else if (ori == "horizontalY")
    {
        return orientationHandling::horizontalY;
    }
    else
    {
        WarningInFunction
            << "Bad well orientation specifier " << ori
            << ", using 'generic'" << endl;
    }
    return orientationHandling::generic;
}


Foam::word Foam::sourceProperties::orientationHandlingToWord
(
    const orientationHandling& ori
)
{
    word enumName;
    switch (ori)
    {
        case orientationHandling::vertical :
        {
            enumName = "vertical";
            break;
        }
        case orientationHandling::horizontalX :
        {
            enumName = "horizontalX";
            break;
        }
        case orientationHandling::horizontalY :
        {
            enumName = "horizontalY";
            break;
        }
        case orientationHandling::generic :
        {
            enumName = "generic";
            break;
        }
    }
    return enumName;
}

Foam::sourceProperties::operationHandling 
Foam::sourceProperties::wordToOpHandling
(
    const word& op
) const
{
    if (op == "production")
    {
        return operationHandling::production;
    }
    else if (op == "injection") 
    {
        return operationHandling::injection;
    }
    else 
    {
        WarningInFunction
            << "Bad well operation mode specifier " << op
            << ", using 'production'" << endl;
    }
    return operationHandling::production;
}


Foam::word Foam::sourceProperties::opHandlingToWord
(
    const operationHandling& op
) const
{
    word enumName("production");
    switch (op)
    {
        case operationHandling::production :
        {
            enumName = "production";
            break;
        }
        case operationHandling::injection :
        {
            enumName = "injection";
            break;
        }
    }
    return enumName;
}

// ************************************************************************* //
