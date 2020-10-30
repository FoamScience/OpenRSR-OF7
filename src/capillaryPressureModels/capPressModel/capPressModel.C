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

#include "regIOobject.H"
#include "capPressModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::string Foam::capPressModel<RockType, nPhases>::nameRegExp = "(.*)<(.*)>";

// * * * * * * * * * * * *  Static Member Functions  * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::Tuple2<Foam::word, Foam::wordList>
Foam::capPressModel<RockType, nPhases>::parseModelName(const word& modelName)
{
    Tuple2<word, wordList> result;
    List<string> groups;

    regExp re(nameRegExp);
    re.match(modelName, groups);

    if (groups.empty())
    {
        WarningInFunction
            << "Can't parse capillary pressure model name: "
            << modelName << nl
            << "Please use simple names like: pcModel<water,oil>" << nl;
        result.first() = modelName;
        return result;
    }
    // Store parsed name
    result.first() = groups[0];
    if (groups.size() == 1)
    {
        return result;
    }

    // Trim spaces
    groups[1].erase
    (
        remove_if(groups[1].begin(), groups[1].end(), ::isspace),
        groups[1].end()
    );

    // Parse comma-separated phases list
    size_t pos = 0;
    while ((pos = groups[1].find(',')) != std::string::npos) {
    result.second().append(groups[1].substr(0, pos));
    groups[1].erase(0, pos + 1);
    }
    result.second().append(groups[1]);

    return result;
}


template<class RockType, int nPhases>
Foam::autoPtr<Foam::capPressModel<RockType, nPhases>>
Foam::capPressModel<RockType, nPhases>::New
(
    const word& name,
    const dictionary& transportProperties,
    const RockType& rock
)
{
    const word modelType =
        transportProperties.subDict(name).lookup("type");

    Info<< "Selecting capillary pressure model type " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown Capillary Pressure Model type"
            << modelType << nl << nl
            << "Valid capillary pressure model types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<capPressModel>
    ( cstrIter()(name, transportProperties, rock) );
}

template<class RockType, int nPhases>
const Foam::capPressModel<RockType, nPhases>&
Foam::capPressModel<RockType, nPhases>::getPcModel
(
    const word& phaseName,
    const fvMesh& mesh
)
{
    word modelName = "";
    HashTable<const capPressModel*> candidates
        = mesh.lookupClass<capPressModel>();
    forAllIter(typename HashTable<const capPressModel*>, candidates, it)
    {
        if (findIndex((*it)->phases(), phaseName) != -1)
        {
            modelName = (*it)->name(); 
        }
    }
    
    if (modelName == "")
    {
        // Make sure the phase exists first
        mesh.lookupObject<phase>(phaseName);
        FatalErrorInFunction
            << "No Pc model for phase " + phaseName + " of type "
            << capPressModel::typeName_() << " was found."
            << exit(FatalError);
    }
    return mesh.lookupObject<capPressModel>(modelName);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::capPressModel<RockType, nPhases>::capPressModel
(
    const word& name,
    const dictionary& transportProperties,
    const RockType& rock
)
:
    regIOobject
    (
        IOobject
        (
            name,
            rock.mesh().time().timeName(),
            rock.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    pcDict_(transportProperties.subDict(name)),
    rock_(rock),
    phaseNames_
    (
        ! parseModelName(name).second().empty()
        ? parseModelName(name).second()
        : pcDict_.lookup("phases")
    ),
    canonicalPhases_
    (
        pcDict_.found("canonicalPhases")
        ? pcDict_.lookup("canonicalPhases")
        : wordList(phaseNames_.begin(),phaseNames_.end()-1)
    ),
    pcTable_(2*canonicalPhases_.size())
{
    if (phaseNames_.size() != nPhases)
    {
        FatalErrorInFunction
            << "Capillary pressure Model " << name
            << " is supposed to have " << nPhases << " phases but got "
            << phaseNames_.size() << nl << nl
            << exit(FatalError);
    }

    // Initialize pcTable
    forAll(canonicalPhases_, pi)
    {
        pcTable_.insert
        (
            pcName(canonicalPhases_[pi]),
            volScalarField
            (
                IOobject
                (
                    pcName(canonicalPhases_[pi]),
                    rock.mesh().time().timeName(),
                    rock.mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                rock.mesh(),
                dimensionedScalar("pc", dimPressure, 0),
                "zeroGradient"
            )
        );
        pcTable_.insert
        (
            dpcName(canonicalPhases_[pi], canonicalPhases_[pi]),
            volScalarField
            (
                IOobject
                (
                    dpcName(canonicalPhases_[pi], canonicalPhases_[pi]),
                    rock.mesh().time().timeName(),
                    rock.mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                rock.mesh(),
                dimensionedScalar("pc", dimPressure, 0),
                "zeroGradient"
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::capPressModel<RockType, nPhases>::~capPressModel()
{}

// ************************************************************************* //
