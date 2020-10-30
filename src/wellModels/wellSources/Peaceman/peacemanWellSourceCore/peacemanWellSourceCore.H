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

Class
    Foam::wellSources::peacemanWellSourceCore

Description
    Number-Of-Phases-Invariant part of the standard peaceman description
    for well source in square grids.

SourceFiles
    peacemanWellSourceCore.C

\*---------------------------------------------------------------------------*/

#ifndef peacemanWellSourceCore_H
#define peacemanWellSourceCore_H

#include "dimensionedScalarFwd.H"
#include "sourceProperties.H"
#include "wellSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace wellSources
{

/*---------------------------------------------------------------------------*\
                    Class peacemanWellSourceCore Declaration
\*---------------------------------------------------------------------------*/

template<class RockType>
class peacemanWellSourceCore
{

protected:

    // Protected data members

        //- Const ref to rock
        const RockType& rock_;

        //- Cell size in the 3 directions
        mutable List<vector> h_;

        //- Cell sizes already estimated ?
        mutable bool estimatedH_;

        //- Equivalent radius
        mutable scalarList re_; 

    // Protected Member functions

        //- Estimate cell sizes
        inline void estimateCellSizes( const labelList& cellIDs );

        //- Estimate requivalent radius
        inline void estimateEquivRadius
        (
            const labelList& cellIDs,
            sourceProperties& srcProps
        );

        //- Calculate Well index
        inline void calculateWellIndex
        (
            const labelList& cellIDs,
            sourceProperties& srcProps
        );


public:

    //- Runtime type information
    ClassName("Peaceman");

    // Constructors

        // Construct from components
        peacemanWellSourceCore
        (
            const word& name,
            const phase& attachedPhase,
            const dictionary& wellSourceDict,
            const RockType& rock
        );

        // Construct from copy
        peacemanWellSourceCore(const peacemanWellSourceCore&) = delete;

	// Destructor
	virtual ~peacemanWellSourceCore();

    // Member operators

        //- Disallow default bitwise assignment
        peacemanWellSourceCore& operator=
        (
            const peacemanWellSourceCore&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wellSources

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "peacemanWellSourceCore.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
