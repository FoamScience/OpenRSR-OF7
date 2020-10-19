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
    Foam::wells::standardWell

Description
    Standard well implementation.

SourceFiles
    standardWell.C

\*---------------------------------------------------------------------------*/

#ifndef standardWell_H
#define standardWell_H

#include "well.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace wells
{

/*---------------------------------------------------------------------------*\
                        Class standardWell Declaration
\*---------------------------------------------------------------------------*/

template<class RockType, int nPhases>
class standardWell
:
    public well<RockType, nPhases>
{
public:

    //- Runtime type information
    ClassName("standard");

    // Constructors

        // Construct from components
        standardWell
        (
            const word& name,
            const dictionary& wellDict,
            const RockType& rock,
            HashTable<autoPtr<wellSource<RockType, nPhases>>>& sources,
            HashTable<fvScalarMatrix>& matTable
        );

        // Construct from copy
        standardWell(const standardWell&) = delete;

	// Destructor
	virtual ~standardWell();

    // Member Functions

        //- Update well sources
        void correct()
        {
            forAll(this->drives_, di)
            {
                this->drives_[di].correct();
            }
        }

    // Member operators

        //- Disallow default bitwise assignment
        standardWell& operator=(const standardWell&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wells

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "standardWell.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
