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

#include "phase.H"
#include "fvcFlux.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcReconstruct.H"
#include "fixedValueFvsPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phase, 0);
    defineRunTimeSelectionTable(phase, dictionary);
}

// * * * * * * * * * * * * Static Members Functions  * * * * * * * * * * * * //

Foam::autoPtr<Foam::phase> Foam::phase::New
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const mixtureType& mT
)
{
    const word phaseType = 
        transportProperties.subDict(name).lookupOrDefault<word>
        (
            "phaseType", "blackoil"
        );

    Info<< "Selecting phase type " << phaseType << " for " << name << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(phaseType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown Phase type "
            << phaseType << nl << nl
            << "Valid Phase types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<phase>
    ( cstrIter()(name, mesh, transportProperties, mT)  );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phase::phase
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const mixtureType& mT
)
:
    objectRegistry
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    phaseDict_(transportProperties.subDict(name)), // Assume it exists
    mesh_(mesh),
    U_
    (
        IOobject
        (
            name+".U",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(name+".U0", dimVelocity, vector::zero)
    ),
    alphaPtr_
    (
        mT == mixtureType::multiPhase
        ?
        new volScalarField
        (
            IOobject
            (
                name+".alpha",
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
        :
        nullptr
    ),
    rhoSc_
    (
        "rhoSc",
        dimDensity, phaseDict_
    ),
    rho0_
    (
        name+".rho0",
        phaseDict_.lookupOrAddDefault<dimensionedScalar>
        (
            "rho0",
            dimensionedScalar("rho0", dimViscosity*dimDensity, 1.0)
        )
    ),
    mu0_
    (
        name+".mu0",
        phaseDict_.lookupOrAddDefault<dimensionedScalar>
        (
            "mu0",
            dimensionedScalar("mu0", dimViscosity*dimDensity, -1)
        )
    ),
    muMesh_
    (
        mu0_.value() != -1
        ?  new singleCellFvMesh
           (
               IOobject
               (
                   name+".muMesh",
                   mesh.polyMesh::instance(),
                   mesh.time(),
                   IOobject::NO_READ,
                   IOobject::NO_WRITE
               ),
               mesh
           )
        :  nullptr
    ),
    phiPtr_
    (
        new surfaceScalarField
        (
           IOobject
           (
               name+".phi",
               mesh.time().timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           linearInterpolate(U_) & mesh.Sf(),
           fixedValueFvsPatchField<vector>::typeName
        )
    ),
    BModel_
    (
        FVFModel::New
        (
            name+".FVFModel", phaseDict_, mesh
        )
    ),
    rho_
    (
        IOobject
        (
            name+".rho",
            mesh.time().timeName(),
            mesh,
            BModel_->isIncompressible()
                ? IOobject::READ_IF_PRESENT
                : IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        BModel_->rFVF().mesh(), // Follow the mesh from BModel
        rho0_
    ),
    mu_
    (
        IOobject
        (
            name+".mu",
            mesh.time().timeName(),
            mesh,
            mu0_.value() == -1
                ? IOobject::MUST_READ
                : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mu0_.value() == -1 ? mesh : muMesh_(),
        dimensionedScalar
        (
            "mu0", dimViscosity*dimDensity,
            mu0_.value() == -1 ? 1.0 : mu0_.value()
        )
    )
{
}


Foam::phase::phase
(
    const phase& ph
)
:
    objectRegistry
    (
        IOobject
        (
            ph.name_,
            ph.mesh_.time().timeName(),
            ph.mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(ph.name_),
    phaseDict_(ph.phaseDict_),
    mesh_(ph.mesh_),
    U_(ph.U_),
    alphaPtr_(ph.alphaPtr_),
    rhoSc_(ph.rhoSc_),
    rho0_(ph.rho0_),
    mu0_(ph.mu0_),
    muMesh_(ph.muMesh_),
    phiPtr_(ph.phiPtr_),
    BModel_(ph.BModel_),
    rho_(ph.rho_),
    mu_(ph.mu_)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phase::~phase() {}

// ************************************************************************* //
