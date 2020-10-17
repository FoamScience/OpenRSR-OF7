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
    Foam::wellModels::peacemanWellModel

Description
    Applies peaceman'ish well source describers for each phase.

SourceFiles
    peacemanWellModel.C

\*---------------------------------------------------------------------------*/

#ifndef peacemanWellModel_H
#define peacemanWellModel_H

#include "wellModel.H"
#include "wellSource.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace wellModels
{

/*---------------------------------------------------------------------------*\
                        Class peacemanWellModel Declaration
\*---------------------------------------------------------------------------*/

template<class RockType>
class peacemanWellModel
:
    public wellModel<RockType, 2>
{
public:

    //- Runtime type information
    ClassName("Peaceman");

    // Constructors

        // Construct from components
        peacemanWellModel
        (
            const word& name,
            const dictionary& transportProperties,
            const dictionary& wellsProperties,
            const RockType& rock
        );

        // Construct from copy
        peacemanWellModel(const peacemanWellModel&) = delete;

	// Destructor
	virtual ~peacemanWellModel();

    // Member Functions

        //- Update well sources
        virtual void correct();

    // Member operators

        //- Disallow default bitwise assignment
        peacemanWellModel& operator=(const peacemanWellModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wellModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "peacemanWellModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
