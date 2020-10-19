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

#include "wellModel.H"
#include "TemplatedRunTimeSelection.H"

namespace Foam {

// Define Types for single-phase well models
defineVariadicTemplatedRunTimeSelectionTable(wellModel, iRock, 1);
// Define Types for two-phase well models
defineVariadicTemplatedRunTimeSelectionTable(wellModel, iRock, 2);
defineVariadicTemplatedRunTimeSelectionTable(wellModel, dRock, 2);
defineVariadicTemplatedRunTimeSelectionTable(wellModel, aRock, 2);

}


// ************************************************************************* //
