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

#include "CFLMethod.H"

// * * * * * * * * * * * *  Static Member Functions  * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::autoPtr<Foam::CFLMethod<RockType, nPhases>>
Foam::CFLMethod<RockType, nPhases>::New
(
    const word& name,
    const dictionary& algorithmProperties,
    const surfaceScalarField& phi,
    const wordList& phaseNames,
    const RockType& rock
)
{
    const word modelType = algorithmProperties.lookup("type");

    Info<< "Selecting CFL Method type " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown CFL Method type"
            << modelType << nl << nl
            << "Valid CFL method types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<CFLMethod>
    ( cstrIter()(name, algorithmProperties,phi, phaseNames, rock) );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::CFLMethod<RockType, nPhases>::CFLMethod
(
    const word& name,
    const dictionary& algorithmProperties,
    const surfaceScalarField& phi,
    const wordList& phaseNames,
    const RockType& rock
)
:
    name_(name),
    cflDict_(algorithmProperties),
    phi_(phi),
    rock_(rock),
    CFLNo_(rock.mesh().nCells(), 0)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::CFLMethod<RockType, nPhases>::~CFLMethod() {}

// ************************************************************************* //
