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

#include "CFLMethod.H"
#include "TemplatedRunTimeSelection.H"

namespace Foam {

// Declare rock types
using iRock = rock<Isotropic>;
using dRock = rock<DiagAnisotropic>;
using aRock = rock<Anisotropic>;

// Define Types for CFL methods
defineVariadicTemplatedRunTimeSelectionTable(CFLMethod, iRock, 2);
defineVariadicTemplatedRunTimeSelectionTable(CFLMethod, dRock, 2);
defineVariadicTemplatedRunTimeSelectionTable(CFLMethod, aRock, 2);

}


// ************************************************************************* //
