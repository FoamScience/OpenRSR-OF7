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

#include "DimensionedField.H"
#include "pcTabular.H"
#include "capPressModel.H"

namespace Foam 
{
namespace twoPhaseCapPressModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType>
pcTabular<RockType>::pcTabular
(
    const word& name,
    const dictionary& transportProperties,
    const RockType& rock
)
:
    capPressModel<RockType, 2>(name, transportProperties, rock),
    alpha_
    (
        rock.mesh().template lookupObject<phase>
            (this->canonicalPhases_[0]).alpha()
    ),
    pcSeries_
    (
        basicInterpolationTable<scalar>::New
        (
            this->pcDict_.subDict("pcData")
        )
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType>
pcTabular<RockType>::~pcTabular() {}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class RockType>
void pcTabular<RockType>::correct()
{
    // Refs to pc fields
    auto& pc = this->pcTable_[this->pcName(this->canonicalPhases_[0])];
    auto& dpc = this->pcTable_
        [this->dpcName(this->canonicalPhases_[0], this->canonicalPhases_[0])];

    forAll(alpha_.internalField(), ci)
    {
        const scalar& alphaInCell = alpha_.internalField()[ci];
        auto interpolatedValues = pcSeries_->interpolate(alphaInCell);
        pc[ci] = interpolatedValues[0];
        dpc[ci] = interpolatedValues[1];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseCapPressModels

} // End namespace Foam

// ************************************************************************* //
