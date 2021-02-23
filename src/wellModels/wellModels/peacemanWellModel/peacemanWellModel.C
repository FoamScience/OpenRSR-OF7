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

#include "error.H"
#include "peacemanWellModel.H"

namespace Foam 
{
namespace wellModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
peacemanWellModel<RockType, nPhases>::peacemanWellModel
(
    const word& name,
    const dictionary& transportProperties,
    const dictionary& wellsProperties,
    const RockType& rock
)
:
    wellModel<RockType,nPhases>(name,transportProperties,wellsProperties,rock)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
peacemanWellModel<RockType, nPhases>::~peacemanWellModel() {}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class RockType, int nPhases>
void peacemanWellModel<RockType, nPhases>::correct()
{
    if (this->wells_.size() == 0)
    {
        return;
    }

    this->clearMatrices();

    forAll(this->wells_, wi)
    {
        this->wells_[wi].correct();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wellModels

} // End namespace Foam

// ************************************************************************* //
