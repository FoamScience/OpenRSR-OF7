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

#include "rock.H"
#include "fixedValueFvsPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class PermeabilityType, class CompressibilityType>
Foam::autoPtr<Foam::rock<PermeabilityType, CompressibilityType>>
Foam::rock<PermeabilityType,CompressibilityType>::New
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& rockProperties
)
{
    // Run with standard rock models by default
    const word rockType = rockProperties.subDict(name).lookupOrDefault<word>
       ("rockType", "standard");

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(rockType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown Rock type "
            << rockType << nl << nl
            << "Valid Rock types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<rock>( cstrIter()(name, mesh, rockProperties) );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class PermeabilityType, class CompressibilityType>
Foam::rock<PermeabilityType, CompressibilityType>::rock
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& rockProperties
)
:
    regIOobject
    (
        IOobject
        (
            name,
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    rockDict_(rockProperties.subDict(name)),
    mesh_(mesh),
    poroInit_("porosity", dimless, rockDict_),
    permInit_("permeability", dimArea, rockDict_),
    compInit_("compressibility", (dimless/dimPressure), rockDict_),
    porosity_
    (
        IOobject
        (
            name+".porosity",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        poroInit_
    ),
    K_
    (
        IOobject
        (
            name+".permeability",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        permInit_
    ),
    c_
    (
        IOobject
        (
            name+".compressibility",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        compInit_
    )
{
}


template<class PermeabilityType, class CompressibilityType>
Foam::rock<PermeabilityType, CompressibilityType>::rock
(
    const rock& rk
)
:
    regIOobject(rk),
    name_(rk.name_),
    rockDict_(rk.rockDict_),
    mesh_(rk.mesh_),
    poroInit_(rk.poroInit_),
    permInit_(rk.permInit_),
    compInit_(rk.compInit_),
    porosity_(rk.porosity_),
    K_(rk.K_),
    c_(rk.c_)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class PermeabilityType, class CompressibilityType>
Foam::rock<PermeabilityType, CompressibilityType>::~rock()
{}

// ************************************************************************* //
