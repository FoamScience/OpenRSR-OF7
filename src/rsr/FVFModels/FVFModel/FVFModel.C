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

#include "FVFModel.H"
#include "dimensionedScalarFwd.H"

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FVFModel, 0);
    defineRunTimeSelectionTable(FVFModel, dictionary);
}

// * * * * * * * * * * * * Static Function Members * * * * * * * * * * * * * //

Foam::autoPtr<Foam::FVFModel> Foam::FVFModel::New
(
    const word& name,
    const dictionary& phaseDict,
    const fvMesh& mesh
)
{
    const word FVFType = phaseDict.lookupOrDefault<word>
    (
        "FVFModel",
        "tabularFVFvsPressure"
    );

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(FVFType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown FVF Model type " << FVFType
            << " for phase " << phaseDict.dictName() << nl << nl
            << "Valid FVF models : " << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<FVFModel>(cstrIter()(name, phaseDict, mesh));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FVFModel::FVFModel
(
    const word& name,
    const dictionary& phaseDict,
    const fvMesh& mesh
)
:
    name_(name),
    phaseDict_(phaseDict),
    mesh_(mesh),
    isIncompressible_
    (
        phaseDict_.lookupOrAddDefault<bool>("incompressible", false)
    ),
    oneCellMesh_
    (
        new singleCellFvMesh
        (
            IOobject
            (
                name+".singleCellMesh",
                mesh.polyMesh::instance(),
                mesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    ),
    rFVF_
    (
        IOobject
        (
            phaseDict.dictName()+".rFVF",
            mesh.time().timeName(),
            mesh,
            isIncompressible_
                ? IOobject::READ_IF_PRESENT
                : IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        isIncompressible_ ? oneCellMesh_() : mesh,
        dimensionedScalar
        (
            "rFVF",
            phaseDict_.lookupOrAddDefault<dimensionedScalar>
            (
                "rFVF", dimensionedScalar("rFVF", dimless, 1.0)
            )
        )
    ),
    drFVFdP_
    (
        IOobject
        (
            phaseDict.dictName()+".drFVFdP",
            mesh.time().timeName(),
            mesh,
            isIncompressible_
                ? IOobject::READ_IF_PRESENT
                : IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        isIncompressible_ ? oneCellMesh_() : mesh,
        dimensionedScalar
        (
            "drFVFdP",
            phaseDict_.lookupOrAddDefault<dimensionedScalar>
            (
                "drFVFdP", dimensionedScalar("drFVFdP", dimless, 0.0)
            )
        )
    )
{}


Foam::FVFModel::FVFModel
(
    const FVFModel& fvfModel
)
:
    name_(fvfModel.name_),
    phaseDict_(fvfModel.phaseDict_),
    mesh_(fvfModel.mesh_),
    isIncompressible_(fvfModel.isIncompressible_),
    oneCellMesh_(fvfModel.oneCellMesh_),
    rFVF_(fvfModel.rFVF_),
    drFVFdP_(fvfModel.drFVFdP_)
{}

// ************************************************************************* //
