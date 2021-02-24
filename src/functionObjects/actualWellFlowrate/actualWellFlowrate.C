/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "actualWellFlowrate.H"
#include "fieldTypes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(actualWellFlowrate, 0);
    addToRunTimeSelectionTable(functionObject, actualWellFlowrate, dictionary);
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::actualWellFlowrate::modeType,
    2
>::names[] = {"reservoir", "surface"};

const Foam::NamedEnum
<
    Foam::functionObjects::actualWellFlowrate::modeType,
    2
> Foam::functionObjects::actualWellFlowrate::modeTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::actualWellFlowrate::writeFileHeader(const label i)
{
    OFstream& file = collated_ ? this->file() : this->files()[i];
    word heading = 
        collated_
        ? word("Actual per-phase Wells flowrate (")
           + modeTypeNames_[mode_] + " conditions)" 
        : names()[i] +
           " actual flowrate (" + modeTypeNames_[mode_] + " conditions)" ;

    writeHeader (file, heading);
    writeCommented(file, "Time");
    collated_ 
        ? writeTabbed(file, "Well\tPhase\tRate")
        : writeTabbed(file, "Rate");

    file<< endl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::actualWellFlowrate::actualWellFlowrate
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    mode_(modeType::reservoir),
    collated_(false),
    wModelName_("wModel"),
    rates_(),
    sourceNotCalculated_(),
    sources_(),
    fileNames_()
{
    read(dict);
    collated_ ? resetName(typeName) : resetNames(fileNames_);
    logFiles::write();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::actualWellFlowrate::~actualWellFlowrate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::actualWellFlowrate::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    mode_ = modeTypeNames_[dict.lookupOrDefault<word>("conditions", "reservoir")];
    collated_ = dict.lookupOrDefault<Switch>("collated", false);
    wModelName_ = dict.lookupOrDefault<word>("wellModel", "wModel");

    const PtrList<entry> wellsInfo ( dict.lookup("wells") );
    forAll(wellsInfo, gi)
    {
        // Select a well group
        const entry& gInfo = wellsInfo[gi];
        if(!gInfo.isDict())
        {
            FatalIOErrorInFunction(dict)
                << "Entry " << gInfo << " in wells list is not a"
                << " valid dictionary." << exit(FatalIOError);
        }
        wordList gPhases (gInfo.dict().lookup("phases"));
        wordList gWellNames (gInfo.dict().lookup("wellNames"));

        // Populate rates table with groups/wells/phases
        forAll(gPhases, pi)
        {
            // Say sources not caclulated yet (at reading stage)
            if (sourceNotCalculated_.found(gPhases[pi]))
            {
                sourceNotCalculated_.set(gPhases[pi], 1);
            } else {
                sourceNotCalculated_.insert(gPhases[pi], 1);
            }

            forAll(gWellNames, wi)
            {
                if (!collated_)
                {
                    // Populate output file names
                    word fName = gWellNames[wi]+"_"+gPhases[pi];
                    if (findIndex(fileNames_, fName) == -1) fileNames_.append(fName);
                }

                // Populate rates Hash table
                rates_.insert
                (
                    gInfo.keyword()+"."+gWellNames[wi]+"."+gPhases[pi],
                    0
                );
            }
        }
    }

    return true;
}


bool Foam::functionObjects::actualWellFlowrate::execute()
{
    return true;
}

bool Foam::functionObjects::actualWellFlowrate::write()
{
    logFiles::write();

    Log << type() << " " << name() <<  " write:" << nl;

    forAll(sourceNotCalculated_.toc(), pi)
    {
        sourceNotCalculated_[sourceNotCalculated_.toc()[pi]] = 1;
    }

    forAll(rates_, wi)
    {
        calcActualRate<iRock, 2>(rates_.toc()[wi], mode_);
        calcActualRate<dRock, 2>(rates_.toc()[wi], mode_);
    }

    sources_.clear();

    forAll(files(), fi) if (Pstream::master()) files()[fi] << endl;
    Log << endl;

    return true;
}



// ************************************************************************* //
