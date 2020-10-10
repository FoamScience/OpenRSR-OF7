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
#include "mathematicalConstants.H"
#include "peacemanWellSource.H"
#include "sourceProperties.H"

namespace Foam 
{
namespace wellSources
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType>
peacemanWellSource<RockType>::peacemanWellSource
(
    const word& name,
    const phase& attachedPhase,
    const dictionary& wellSourceDict,
    const RockType& rock
)
:
    wellSource<RockType, 2>(name, attachedPhase, wellSourceDict, rock),
    h_(),
    estimatedH_(false),
    re_()
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType>
peacemanWellSource<RockType>::~peacemanWellSource() {}

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class RockType>
void peacemanWellSource<RockType>::estimateCellSizes(const labelList& cellIDs)
{
    h_.resize(cellIDs.size());
    // First estimate cell sizes
    forAll(h_, ci)
    {
        const label& cellID = cellIDs[ci];
        const fvMesh& mesh = this->rock_.mesh();
        const labelList& cellEdges = mesh.cellEdges()[cellID];

        scalarList dEdgesAlongX(cellEdges.size());
        scalarList dEdgesAlongY(cellEdges.size());
        scalarList dEdgesAlongZ(cellEdges.size());

        forAll(cellEdges, edgei){
            dEdgesAlongX[edgei] = 
                mag( mesh.edges()[cellEdges[edgei]].vec(mesh.points()).x() );
            dEdgesAlongY[edgei] = 
                mag( mesh.edges()[cellEdges[edgei]].vec(mesh.points()).y() );
            dEdgesAlongZ[edgei] = 
                mag( mesh.edges()[cellEdges[edgei]].vec(mesh.points()).z() );
        }
        h_[ci].x() = gMax(dEdgesAlongX);
        h_[ci].y() = gMax(dEdgesAlongY);
        h_[ci].z() = gMax(dEdgesAlongZ);
    }
    estimatedH_ = true;
}

template<class RockType>
void peacemanWellSource<RockType>::estimateEquivRadius
(
    const labelList& cellIDs,
    sourceProperties& srcProps
)
{
    // Equivalent radius depends on Rock Type and well direction;
    // So, do not implement the generale case and rely on partial
    // specializations
    NotImplemented;
}

template<class RockType>
void peacemanWellSource<RockType>::calculateWellIndex
(
    const labelList& cellIDs,
    sourceProperties& srcProps
)
{
    // Well index depends on equivalent radius
    NotImplemented;
}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class RockType>
void peacemanWellSource<RockType>::calculateCoeff0
(
    scalarList &coeff0,
    sourceProperties& srcProps,
    const labelList& cellIDs
)
{
    word krName = relPermModel<RockType,2>::krName(this->phase_.name());
    const auto& kr = this->krModel_[krName];
    const auto& mu = this->phase_.mu();

    calculateWellIndex(cellIDs, srcProps);
    coeff0.resize(srcProps.wellIndex().size());
    forAll(coeff0, ci)
    {
        const label cellID = cellIDs[ci];
        coeff0[ci] = - srcProps.wellIndex()[ci] * kr[cellID] / mu[cellID];
    }
}

template<class RockType>
void peacemanWellSource<RockType>::calculateCoeff1
(
    scalarList &coeff1,
    sourceProperties& srcProps,
    const labelList& cellIDs
)
{
    word krName = relPermModel<RockType,2>::krName(this->phase_.name());
    const auto& kr = this->krModel_[krName];
    const auto& mu = this->phase_.mu();

    calculateWellIndex(cellIDs, srcProps);
    coeff1.resize(srcProps.wellIndex().size());
    forAll(coeff1, ci)
    {
        const label cellID = cellIDs[ci];
        coeff1[ci] = srcProps.wellIndex()[ci] * kr[cellID] / mu[cellID];
    }
}

template<class RockType>
void peacemanWellSource<RockType>::calculateCoeff2
(
    scalarList &coeff2,
    sourceProperties& srcProps,
    const labelList& cellIDs
)
{
    //const auto& rho = this->phase_.rho();

    calculateWellIndex(cellIDs, srcProps);
    coeff2.resize(srcProps.wellIndex().size());
    // TODO: Add gravity support
    forAll(coeff2, ci)
    {
        coeff2[ci] = 0;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wellSources

} // End namespace Foam

#include "peacemanWellSourceIsoSpec.C"

// ************************************************************************* //
