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
#include "regIOobject.H"
#include "relPermModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::string Foam::relPermModel<RockType, nPhases>::nameRegExp = "(.*)<(.*)>";

// * * * * * * * * * * * *  Static Member Functions  * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::Tuple2<Foam::word, Foam::wordList>
Foam::relPermModel<RockType, nPhases>::parseModelName(const word& modelName)
{
    Tuple2<word, wordList> result;
    List<string> groups;

    regExp re(nameRegExp);
    re.match(modelName, groups);

    if (groups.empty())
    {
        WarningInFunction
            << "Can't parse relative permeability model name: "
            << modelName << nl
            << "Please use simple names like: krModel<water, oil>" << nl;
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
Foam::autoPtr<Foam::relPermModel<RockType, nPhases>>
Foam::relPermModel<RockType, nPhases>::New
(
    const word& name,
    const dictionary& transportProperties,
    const RockType& rock
)
{
    const word modelType =
        transportProperties.subDict(name).lookup("type");

    Info<< "Selecting relative permeability model type " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown Relative Permeability Model type"
            << modelType << nl << nl
            << "Valid relative permeability model types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<relPermModel>
    ( cstrIter()(name, transportProperties, rock) );
}


template<class RockType, int nPhases>
const Foam::relPermModel<RockType, nPhases>&
Foam::relPermModel<RockType, nPhases>::getKrModel
(
    const word& phaseName,
    const fvMesh& mesh
)
{
    word modelName = "";
    HashTable<const relPermModel*> candidates
        = mesh.lookupClass<relPermModel>();
    forAllIter(typename HashTable<const relPermModel*>, candidates, it)
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
            << "No Kr model for phase " + phaseName + " of type "
            << relPermModel::typeName_() << " was found."
            << exit(FatalError);
    }
    return mesh.lookupObject<relPermModel>(modelName);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::relPermModel<RockType, nPhases>::relPermModel
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
    krDict_(transportProperties.subDict(name)),
    rock_(rock),
    phaseNames_
    (
        ! parseModelName(name).second().empty()
        ? parseModelName(name).second()
        : krDict_.lookup("phases")
    ),
    canonicalPhases_
    (
        krDict_.found("canonicalPhases")
        ? krDict_.lookup("canonicalPhases")
        : wordList(phaseNames_.begin(),phaseNames_.end()-1)
    ),
    krTable_(nPhases*(1+canonicalPhases_.size()))
{
    if (phaseNames_.size() != nPhases)
    {
        FatalErrorInFunction
            << "Relative permeability Model " << name
            << " is supposed to have " << nPhases << " phases but got "
            << phaseNames_.size() << nl << nl
            << exit(FatalError);
    }

    // Initialize krTable
    // Reserving space for a kr field per phase, and then for derivatives
    // wrt canonical phases per phase
    forAll(phaseNames_, pi)
    {
        krTable_.insert
        (
            krName(phaseNames_[pi]),
            new volScalarField
            (
                IOobject
                (
                    krName(phaseNames_[pi]),
                    rock.mesh().time().timeName(),
                    rock.mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                rock.mesh(),
                dimensionedScalar("kr", dimless, 1),
                "zeroGradient"
            )
        );
        forAll(canonicalPhases_, ci)
        {
            krTable_.insert
            (
                dkrName(phaseNames_[pi], canonicalPhases_[ci]),
                new volScalarField
                (
                    IOobject
                    (
                        dkrName(phaseNames_[pi], canonicalPhases_[ci]),
                        rock.mesh().time().timeName(),
                        rock.mesh(),
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    rock.mesh(),
                    dimensionedScalar("kr", dimless, 0),
                    "zeroGradient"
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::relPermModel<RockType, nPhases>::~relPermModel()
{}

// ************************************************************************* //
