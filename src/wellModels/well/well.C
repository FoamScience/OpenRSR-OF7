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

#include "UList.H"
#include "error.H"
#include "regIOobject.H"
#include "well.H"

// * * * * * * * * * * * *  Static Member Functions  * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::autoPtr<Foam::well<RockType, nPhases>>
Foam::well<RockType, nPhases>::New
(
    const word& name,
    const dictionary& wellDict,
    const RockType& rock
)
{
    const word modelType = wellDict.lookup("type");

    Info<< tab << "Selecting well type " << modelType << " for "
        << name << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown Well type"
            << modelType << nl << nl
            << "Valid well types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<well>
    ( cstrIter()(name, wellDict, rock) );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::well<RockType, nPhases>::well
(
    const word& name,
    const dictionary& wellDict,
    const RockType& rock
)
:
    List<autoPtr<regIOobject> >(),
    name_(wellDict.dictName()),
    wellDict_(wellDict),
    rock_(rock),
    wellGroups_
    (
        wellDict.lookupOrDefault<wordList>("groups", wordList(1, "defaultGrp"))
    ),
    operation_
    (
        wordToOpHandling(wellDict.lookup("operationMode"))
    ),
    srcProps_(rock.mesh(), wellDict),
    injPhase_
    (
     operation_ == operationHandling::injection
     ? wellDict.lookupOrDefault<word>("injectedPhase", "water")
     : "none"
    ),
    perfos_()
{
    registerToGroups();
    readPerforations();
    readImposedDrives();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::well<RockType, nPhases>::~well()
{}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class RockType, int nPhases>
typename Foam::well<RockType, nPhases>::operationHandling 
Foam::well<RockType, nPhases>::wordToOpHandling
(
    const word& op
) const
{
    if (op == "production")
    {
        return operationHandling::production;
    }
    else if (op == "injection") 
    {
        return operationHandling::injection;
    }
    else 
    {
        WarningInFunction
            << "Bad well operation mode specifier " << op
            << ", using 'production'" << endl;
    }
    return operationHandling::production;
}


template<class RockType, int nPhases>
Foam::word Foam::well<RockType, nPhases>::opHandlingToWord
(
    const operationHandling& op
) const
{
    word enumName("production");
    switch (op)
    {
        case operationHandling::production :
        {
            enumName = "production";
            break;
        }
        case operationHandling::injection :
        {
            enumName = "injection";
            break;
        }
    }
    return enumName;
}

template<class RockType, int nPhases>
void Foam::well<RockType, nPhases>::registerToGroups()
{
    this->List<autoPtr<regIOobject>>::setSize(wellGroups_.size());
    forAll(wellGroups_, gi)
    {
        const objectRegistry& wg =
            rock_.mesh().template lookupObject<objectRegistry>
            (
                wellGroups_[gi]
            );
        this->List<autoPtr<regIOobject>>::operator[](gi).set
        (
            new wellGrpIO
            (
                IOobject
                (
                    wellGroups_[gi]+name_,
                    rock_.mesh().time().timeName(),
                    wg,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}

template<class RockType, int nPhases>
void Foam::well<RockType, nPhases>::readPerforations()
{
    Info << tab << "Constructing cells for well: " << name_ << nl;

    // Perforated intervals as entries
    const PtrList<entry> perfsInfo ( wellDict_.lookup("perforations") );
    // Reshape perforations list
    perfos_.setSize(perfsInfo.size());

    // Construct cells from given topoSetSources
    forAll(perfos_, perfi)
    {
        // Select a perforation interval
        const entry& perfInfo = perfsInfo[perfi];
        // Require that the perforation interval is a valid dict
        if(!perfInfo.isDict())
        {
            FatalIOErrorInFunction(wellDict_)
                << "Entry " << perfInfo << " in wells section is not a"
                << " valid dictionary." << exit(FatalIOError);
        }
        // Set the pointer to the requested topoSetSource
        perfos_.set(
            perfi,
            topoSetSource::New(perfInfo.keyword(),rock_.mesh(),perfInfo.dict())
        );
        // Include selected cells in the well's cell set
        perfos_[perfi].applyToSet(topoSetSource::ADD, wellSet_);
    }

    // Write well set to disk
    if (debug) wellSet_.write();
}


template<class RockType, int nPhases>
void Foam::well<RockType, nPhases>::readImposedDrives()
{
    Info << tab << "Constructing drives for well: " << name_ << nl;

    // Drives as entries
    const PtrList<entry> drivesInfo ( wellDict_.lookup("imposedDrives") );
    // Reshape perforations list
    drives_.setSize(drivesInfo.size());

    // Construct drive objects from cells
    forAll(drives_, di)
    {
        // Select a drive
        const entry& driveInfo = drivesInfo[di];
        // Require that the drive entry is a valid dict
        if(!driveInfo.isDict())
        {
            FatalIOErrorInFunction(wellDict_)
                << "Entry " << driveInfo << " in wells section is not a"
                << " valid dictionary." << exit(FatalIOError);
        }
        // Set the pointer to the requested topoSetSource
        // TODO: Probably use driveHandling::New
        //perfos_.set
        //(
        //    di,
        //    topoSetSource::New
        //    (
        //        driveInfo.keyword(),
        //        rock_.mesh(),
        //        driveInfo.dict()
        //    )
        //);
    }
}
// ************************************************************************* //
