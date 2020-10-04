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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sourceProperties::sourceProperties
(
    const fvMesh& mesh,
    const dictionary& wellDict
)
:
    g_
    (
        mesh.lookupObject<UniformDimensionedField<scalar>>("g")
    ),
    orientation_
    (
        wordToOrientationHandling
        (
            wellDict.lookupOrDefault<word>("orientation", "vertical")
        )
    ),
    radius_
    (
        wellDict.dictName()+".radius",
        dimensionedScalar("radius", dimLength, wellDict)
    ),
    skin_
    (
        readScalar(wellDict.lookup("skin"))
    ),
    J_()
{
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
    else 
    {
        WarningInFunction
            << "Bad well orientation mode specifier " << ori
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
        case orientationHandling::generic :
        {
            enumName = "generic";
            break;
        }
    }
    return enumName;
}

// ************************************************************************* //
