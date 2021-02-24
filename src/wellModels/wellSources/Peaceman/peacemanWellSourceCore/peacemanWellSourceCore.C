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

#include "mathematicalConstants.H"
#include "peacemanWellSourceCore.H"
#include "sourceProperties.H"

namespace Foam 
{
namespace wellSources
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType>
peacemanWellSourceCore<RockType>::peacemanWellSourceCore
(
    const word& name,
    const phase& attachedPhase,
    const dictionary& wellSourceDict,
    const RockType& rock
)
:
    rock_(rock),
    h_(),
    estimatedH_(false),
    re_()
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType>
peacemanWellSourceCore<RockType>::~peacemanWellSourceCore() {}

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class RockType>
inline void peacemanWellSourceCore<RockType>::estimateCellSizes
(
    const labelList& cellIDs
)
{
    h_.resize(cellIDs.size());
    // First estimate cell sizes
    forAll(h_, ci)
    {
        const label& cellID = cellIDs[ci];
        const fvMesh& mesh = rock_.mesh();
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
inline void peacemanWellSourceCore<RockType>::estimateEquivRadius
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
inline void peacemanWellSourceCore<RockType>::calculateWellIndex
(
    const labelList& cellIDs,
    sourceProperties& srcProps
)
{
    // Well index depends on equivalent radius
    NotImplemented;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wellSources

} // End namespace Foam

#ifdef NoRepository
    #include "peacemanWellSourceCoreIsoSpec.C"
#endif

// ************************************************************************* //
