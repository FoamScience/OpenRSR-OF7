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

#include "ListOps.H"
#include "wellModel.H"
#include <set>

// * * * * * * * * * * * *  Static Member Functions  * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::autoPtr<Foam::wellModel<RockType, nPhases>>
Foam::wellModel<RockType, nPhases>::New
(
    const word& name,
    const dictionary& wellsProperties,
    const RockType& rock
)
{
    const word modelType = wellsProperties.lookup("wellModel");

    Info<< "Selecting well model type " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown Well Model Model type"
            << modelType << nl << nl
            << "Valid well model types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<wellModel>
    ( cstrIter()(name, wellsProperties, rock) );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::wellModel<RockType, nPhases>::wellModel
(
    const word& name,
    const dictionary& wellsProperties,
    const RockType& rock
)
:
    name_(name),
    wellsProperties_(wellsProperties),
    rock_(rock),
    groups_(),
    wells_(),
    sources_(nPhases)
{
    createWells();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::wellModel<RockType, nPhases>::~wellModel()
{}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::PtrList<Foam::entry>
Foam::wellModel<RockType, nPhases>::assembleWellGroups()
{
    std::set<word> grpNames;
    // Wells entries in dict
    const PtrList<entry> wells( wellsProperties_.lookup("wells") );
    // Go through each entry and collect mentioned groups
    forAll(wells, wi)
    {
        const entry& we = wells[wi];
        if(!we.isDict())
        {
            FatalIOErrorInFunction(wellsProperties_)
                << "Entry " << we << "in wells section is not a valid "
                << "dictionary." << exit(FatalIOError);
        }
        wordList names = we.dict().lookupOrDefault<wordList>
        (
            "groups", wordList(1, "defaultGrp")
        );
        forAll(names, ni)
        {
            grpNames.insert(names[ni]);
        }
    }
    groups_.setSize(grpNames.size());
    int gi = 0;
    forAllConstIter(std::set<word>, grpNames, it)
    {
        groups_.set
        (
            gi,
            new objectRegistry
            (
                IOobject
                (
                    *it,
                    rock_.mesh().time().constant(),
                    rock_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                )
            )
        );
        gi++;
    }
    return wells;
}

template<class RockType, int nPhases>
void Foam::wellModel<RockType, nPhases>::createWells()
{
    // Parse well entries and create well groups 
    PtrList<entry> wellEntries = assembleWellGroups();
    wells_.setSize(wellEntries.size());

    // Construct well objects
    forAll(wells_, wi)
    {
        const dictionary& wDict = wellEntries[wi].dict();
        wells_.set
        (
            wi,
            well<RockType, nPhases>::New(wDict.dictName(), wDict, rock_)
        );
    }
}

// ************************************************************************* //
