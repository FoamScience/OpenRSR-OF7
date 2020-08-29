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
#include "fixedValueFvsPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompressibilityType, class ViscosityType>
Foam::phase<CompressibilityType, ViscosityType>::phase
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const mixtureType& mT
)
:
    regIOobject
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
    phaseDict_(dict.subDict(name)),
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
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(name+".alpha", dimless, 1.0)
        )
        :
        nullptr
    ),
    rhoSc_
    (
        name+".rhoSc",
        dimDensity,
        readScalar(phaseDict_.lookup("rhoSc"))
    ),
    rho_
    (
        IOobject
        (
            name+".rho",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        rhoSc_
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
    mu_
    (
        IOobject
        (
            name+".rho",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(name+".rho0", dimless, 1.0) 
    ),
    BModel_
    (
        FVFModel<CompressibilityType>::New
        (
            name+".FVFModel", phaseDict_, mesh
        )
    )
{
}


template<class CompressibilityType, class ViscosityType>
Foam::phase<CompressibilityType, ViscosityType>::phase
(
    const phase& ph
)
:
    regIOobject(ph),
    name_(ph.name_),
    phaseDict_(ph.phaseDict_),
    mesh_(ph.mesh_),
    U_(ph.U_),
    alphaPtr_(ph.alphaPtr_),
    rho_(ph.rho_),
    phiPtr_(ph.phiPtr_),
    mu_(ph.mu_),
    BModel_
    (
        FVFModel<CompressibilityType>::New
        (
            name_+".FVFModel", phaseDict_, mesh_
        )
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompressibilityType, class ViscosityType>
Foam::phase<CompressibilityType, ViscosityType>::~phase()
{}

// ************************************************************************* //
