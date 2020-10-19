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

#include "dimensionedScalarFwd.H"
#include "rock.H"
#include "fixedValueFvsPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class PermeabilityType>
Foam::autoPtr<Foam::rock<PermeabilityType>>
Foam::rock<PermeabilityType>::New
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

template<class PermeabilityType>
Foam::rock<PermeabilityType>::rock
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
    isIncompressible_(rockDict_.lookupOrAddDefault("incompressible", false)),
    poroInit_
    (
        rockDict_.found("porosity")
        ? dimensionedScalar("porosity", dimless, rockDict_)
        : dimensionedScalar("porosity", dimless, -1)
    ),
    compInit_
    (
        rockDict_.found("compressibility")
        ? dimensionedScalar("compressibility", dimless/dimPressure, rockDict_)
        : dimensionedScalar("compressibility", dimless/dimPressure, -1)
    ),
    permInit_
    (
        rockDict_.found("permeability")
        ? dimensioned<KcmptType>("permeability", dimArea, rockDict_)
        : dimensioned<KcmptType>
          ("permeability", dimArea, static_cast<KcmptType>(Zero))
    ),
    oneCellMesh_
    (
        ( compInit_.value() == 0 and poroInit_.value() == 0 )
        ? nullptr
        : new singleCellFvMesh
        (
            IOobject
            (
                name+".oneCellMesh",
                mesh.polyMesh::instance(),
                mesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    ),
    porosity_
    (
        poroInit_.value() == -1
        ?   // If no single value is provided in dict
            std::move
            (
             VolatileDimensionedField<scalar,volMesh>
                (
                    IOobject
                    (
                        name+".porosity",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            )
        :   // If a single value is provided
            std::move
            (
             VolatileDimensionedField<scalar,volMesh>
                (
                    IOobject
                    (
                        name+".porosity",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    oneCellMesh_,
                    poroInit_
                )
            )
    ),
    c_
    (
        compInit_.value() == -1
        ?   // If no single value is provided in dict
            std::move
            (
             VolatileDimensionedField<scalar,volMesh>
                (
                    IOobject
                    (
                        name+".compressibility",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            )
        :   // If a single value is provided
            std::move
            (
             VolatileDimensionedField<scalar,volMesh>
                (
                    IOobject
                    (
                        name+".compressibility",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    oneCellMesh_,
                    compInit_
                )
            )
    ),
    K_
    (
        permInit_.value() == static_cast<KcmptType>(Zero)
        ?
            PermeabilityType
            (
                IOobject
                (
                    name+".permeability",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        :
            PermeabilityType
            (
                IOobject
                (
                    name+".permeability",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                permInit_,
                "zeroGradient"
            )
    )
{
}


template<class PermeabilityType>
Foam::rock<PermeabilityType>::rock
(
    const rock& rk
)
:
    regIOobject(rk),
    name_(rk.name_),
    rockDict_(rk.rockDict_),
    mesh_(rk.mesh_),
    isIncompressible_(rk.isIncompressible_),
    poroInit_(rk.poroInit_),
    compInit_(rk.compInit_),
    permInit_(rk.permInit_),
    oneCellMesh_(rk.oneCellMesh_),
    porosity_(rk.porosity_),
    c_(rk.c_),
    K_(rk.K_)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class PermeabilityType>
Foam::rock<PermeabilityType>::~rock()
{}

// ************************************************************************* //
