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
    Foam::wellSources::peacemanWellSource

Description
    Peaceman's description for the relationship q = f(p) in a square grid:
    q_ = J*kr/mu(BHP - p_ - gamma (BHZ - z))

SourceFiles
    peacemanWellSource.C

\*---------------------------------------------------------------------------*/

#ifndef peacemanWellSource_H
#define peacemanWellSource_H

#include "dimensionedScalarFwd.H"
#include "sourceProperties.H"
#include "wellSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace wellSources
{

/*---------------------------------------------------------------------------*\
                        Class peacemanWellSource Declaration
\*---------------------------------------------------------------------------*/

template<class RockType>
class peacemanWellSource
:
    public wellSource<RockType, 2>
{

protected:

    // Protected data members

        //- Cell size in the 3 directions
        mutable List<vector> h_;

        //- Cell sizes already estimated ?
        mutable bool estimatedH_;

        //- Equivalent radius
        mutable scalarList re_; 

    // Protected Member functions

        //- Estimate cell sizes
        void estimateCellSizes( const labelList& cellIDs );

        //- Estimate requivalent radius
        void estimateEquivRadius
        (
            const labelList& cellIDs,
            sourceProperties& srcProps
        );

        //- Calculate Well index
        void calculateWellIndex
        (
            const labelList& cellIDs,
            sourceProperties& srcProps
        );


public:

    //- Runtime type information
    ClassName("Peaceman");

    // Constructors

        // Construct from components
        peacemanWellSource
        (
            const word& name,
            const phase& attachedPhase,
            const dictionary& wellSourceDict,
            const RockType& rock
        );

        // Construct from copy
        peacemanWellSource(const peacemanWellSource&) = delete;

	// Destructor
	virtual ~peacemanWellSource();

    // Member Functions

        //- How to calculate coeff0 in q = coeff0 p + coeff1 BHP + coeff2
        virtual void calculateCoeff0
        (
            scalarList& coeff0,
            sourceProperties& srcProps,
            const labelList& cellIDs
        );

        //- How to calculate coeff1 in q = coeff0 p + coeff1 BHP + coeff2
        virtual void calculateCoeff1
        (
            scalarList& coeff1,
            sourceProperties& srcProps,
            const labelList& cellIDs
        );

        //- How to calculate coeff2 in q = coeff0 p + coeff1 BHP + coeff2
        virtual void calculateCoeff2
        (
            scalarList& coeff2,
            sourceProperties& srcProps,
            const labelList& cellIDs
        );


    // Member operators

        //- Disallow default bitwise assignment
        peacemanWellSource& operator=(const peacemanWellSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wellSources

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "peacemanWellSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
