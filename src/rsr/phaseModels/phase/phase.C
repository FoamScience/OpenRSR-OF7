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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CompressibilityType, class ViscosityType>
Foam::autoPtr<Foam::phase<CompressibilityType, ViscosityType>>
Foam::phase<CompressibilityType, ViscosityType>::New
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const mixtureType& mT
)
{
    const word phaseType = 
        transportProperties.subDict(name).lookup("phaseType");

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

template<class CompressibilityType, class ViscosityType>
Foam::phase<CompressibilityType, ViscosityType>::phase
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const mixtureType& mT
)
:
    phaseCore(name, mesh, transportProperties, mT),
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
    phaseCore(ph),
    rho_(ph.rho_),
    mu_(ph.mu_),
    BModel_ // Creates new model
    (
        FVFModel<CompressibilityType>::New
        (
            this->name_+".FVFModel", this->phaseDict_, this->mesh_
        )
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompressibilityType, class ViscosityType>
Foam::phase<CompressibilityType, ViscosityType>::~phase()
{}

// ************************************************************************* //
