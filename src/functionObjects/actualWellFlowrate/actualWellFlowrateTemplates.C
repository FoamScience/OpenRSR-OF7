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
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
void Foam::functionObjects::actualWellFlowrate::calcActualRate
(
    const word& rateName,
    const modeType& mode
)
{
    // TODO: Parallelize this code
    word nameCopy = rateName;

    // Find the well model, else exit
    if (!mesh_.foundObject<wellModel<RockType, nPhases>>(wModelName_))
    {
        if (debug)
        {
            const wellModel<RockType, nPhases>& wModel =
                mesh_.lookupObject<wellModel<RockType, nPhases>>(wModelName_);
        }
        return;
    }
    const wellModel<RockType, nPhases>& wModel =
        mesh_.lookupObject<wellModel<RockType, nPhases>>(wModelName_);

    // Split well name
    int dpos = nameCopy.find('.');
    word gName = nameCopy.substr(0, dpos);
    nameCopy.erase(0, dpos+1);
    dpos = nameCopy.find('.');
    word wName = nameCopy.substr(0, dpos);
    nameCopy.erase(0, dpos+1);
    dpos = nameCopy.find('.');
    word pName = nameCopy.substr(0, dpos);
    nameCopy.erase(0, dpos+1);

    // Prepare file to write to
    label fIndex = findIndex(names(), wName+"_"+pName);
    OFstream& file = collated_ ? this->file() : this->files()[fIndex];

    // Find the group/well
    const objectRegistry& wg = mesh_.lookupObject<objectRegistry>(gName);
    const well<RockType, nPhases>& w =
        wg.lookupObject<well<RockType, nPhases>>(wName);

    // Calculate the rate
    if (sourceNotCalculated_[pName])
    {
        sources_.insert(pName, wModel.explicitSource(pName));
        sourceNotCalculated_[pName] = 0;
    }

    scalar rate = 0;
    forAll(w.cellIDs(), ci)
    {
        rate += sources_[pName][ci];
    }

    writeTime(file);
    if (collated_)
    {
        file<< tab << wName << tab << pName << tab << mag(rate) << endl; 
    } else {
        file<< tab << mag(rate);
    }

    Log<< gName << " group, " << wName << " well, " << pName << " rate : "
       << rate << endl;

    rates_.set(rateName, rate);

}


// ************************************************************************* //
