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
#include "krTabular.H"
#include "oneField.H"
#include "relPermModel.H"

namespace Foam 
{
namespace twoPhaseRelPermModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType>
krTabular<RockType>::krTabular
(
    const word& name,
    const dictionary& transportProperties,
    const RockType& rock
)
:
    relPermModel<RockType, 2>(name, transportProperties, rock),
    otherPhase_
    (
        this->phaseNames_[0] == this->canonicalPhases_[0]
        ? this->phaseNames_[1]
        : this->phaseNames_[0]
    ),
    alpha_
    (
        rock.mesh().template lookupObject<phase>
            (this->canonicalPhases_[0]).alpha()
    ),
    krSeries_
    (
        basicInterpolationTable<scalar>::New
        (
            this->krDict_.subDict("krData")
        )
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType>
krTabular<RockType>::~krTabular() {}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class RockType>
void krTabular<RockType>::correct()
{
    // Refs to kr fields
    auto& kr1 = this->operator[](this->krName(this->canonicalPhases_[0]));
    auto& kr2 = this->operator[](this->krName(otherPhase_));
    auto& dkr1 = this->operator[]
        (this->dkrName(this->canonicalPhases_[0], this->canonicalPhases_[0]));
    auto& dkr2 = this->operator[]
        (this->dkrName(otherPhase_, this->canonicalPhases_[0]));

    forAll(alpha_.internalField(), ci)
    {
        const scalar& alphaInCell = alpha_.internalField()[ci];
        auto interpolatedValues = krSeries_->interpolate(alphaInCell);
        kr1[ci] = interpolatedValues[0];
        kr2[ci] = interpolatedValues[1];
        dkr1[ci] = interpolatedValues[2];
        dkr2[ci] = interpolatedValues[3];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseRelPermModels

} // End namespace Foam

// ************************************************************************* //
