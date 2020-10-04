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

#include "volFieldsFwd.H"
#include "wellSource.H"

// * * * * * * * * * * * *  Static Member Functions  * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::autoPtr<Foam::wellSource<RockType, nPhases>>
Foam::wellSource<RockType, nPhases>::New
(
    const word& name,
    const phase& attachedPhase,
    const dictionary& wellSourceDict,
    const RockType& rock
)
{
    const word modelType = wellSourceDict.lookup("wellSourceType");

    Info<< "Selecting well source type " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown Well Source type"
            << modelType << nl << nl
            << "Valid well source types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<wellSource>
    ( cstrIter()(name, attachedPhase, wellSourceDict, rock) );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::wellSource<RockType, nPhases>::wellSource
(
    const word& name,
    const phase& attachedPhase,
    const dictionary& wellSourceDict,
    const RockType& rock
)
:
    name_(name),
    phase_(attachedPhase),
    wellSourceDict_(wellSourceDict),
    rock_(rock),
    p_(rock.mesh().template lookupObject<volScalarField>("p")),
    krModel_
    (
        relPermModel<RockType, nPhases>::getKrModel
        (
            attachedPhase.name(), rock.mesh()
        )
    ),
    pcModel_
    (
        capPressModel<RockType, nPhases>::getPcModel
        (
            attachedPhase.name(), rock.mesh()
        )
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
Foam::wellSource<RockType, nPhases>::~wellSource()
{}

// ************************************************************************* //
