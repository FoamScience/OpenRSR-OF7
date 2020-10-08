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

#include "peacemanWellSource.H"
#include "sourceProperties.H"

namespace Foam 
{
namespace wellSources
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<>
void peacemanWellSource<iRock>::estimateEquivRadius
(
    const labelList& cellIDs,
    sourceProperties& srcProps
)
{
    // Estimate cell sizes if not already done
    if (!estimatedH_) estimateCellSizes(cellIDs);
    re_.resize(cellIDs.size());

    // TODO: Fix orientation handling
    direction id1 = 
        srcProps.orientation() 
        == sourceProperties::orientationHandling::generic ? 1 : 0;
    direction id2 = 
        srcProps.orientation() 
        == sourceProperties::orientationHandling::vertical ? 1 : 2;
        
    // Isotropic medium --> re depends only on geometry
    forAll(re_, ci)
    {
        // TODO: Adjust to well direction
        re_[ci] = 
            0.14*
            Foam::sqrt(Foam::pow(h_[ci][id1],2) + Foam::pow(h_[ci][id2],2));
    }
}

template<>
void peacemanWellSource<iRock>::calculateWellIndex
(
    const labelList cellIDs,
    sourceProperties& srcProps
)
{
    estimateEquivRadius(cellIDs, srcProps);
    scalarList& J = srcProps.wellIndex();
    J.resize(cellIDs.size());
    auto id = 
        srcProps.orientation() 
        == sourceProperties::orientationHandling::vertical ? 2 : 0;

    forAll(J, ci)
    {
        label cellID = cellIDs[ci];
        J[ci] =
           2 * constant::mathematical::pi * this->rock_.K()[cellID] 
           * h_[ci][id]
           / (log(re_[ci]/srcProps.radius().value()) + srcProps.skin());
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wellSources

} // End namespace Foam

// ************************************************************************* //
