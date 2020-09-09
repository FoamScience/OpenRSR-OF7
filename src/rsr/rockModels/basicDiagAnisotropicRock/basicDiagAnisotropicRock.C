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
    along with OpenFOAM.  If not, see .

\*---------------------------------------------------------------------------*/

#include "basicDiagAnisotropicRock.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace rocks {
        defineTypeNameAndDebug (basicDiagAnisotropicRock, 0);
        addToTemplateRockRunTimeSelectionTable
        (
            rock, basicDiagAnisotropicRock, DiagAnisotropic, dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rocks::basicDiagAnisotropicRock::
basicDiagAnisotropicRock
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& rockProperties
)
:
    rock(name, mesh, rockProperties)
{
}


Foam::rocks::basicDiagAnisotropicRock::
basicDiagAnisotropicRock
(
    const basicDiagAnisotropicRock& rk
)
:
    rock(rk)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rocks::basicDiagAnisotropicRock::~basicDiagAnisotropicRock() {}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

void Foam::rocks::basicDiagAnisotropicRock::correct()
{
    // TODO: Do absolutely nothing for now
    // TODO: Or maybe update compressibility based on pressure
    // Alternative: Update compressibility using Hall's Correlation
}

// ************************************************************************* //
