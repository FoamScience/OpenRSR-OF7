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

// * * * * * * * * * * * * Static Function Members * * * * * * * * * * * * * //

template<class CompressibilityType>
Foam::autoPtr<Foam::FVFModel<CompressibilityType>>
Foam::FVFModel<CompressibilityType>::New
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
            << nl << nl
            << "Valid FVF models : " << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<FVFModel>(cstrIter()(name, phaseDict, mesh));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompressibilityType>
Foam::FVFModel<CompressibilityType>::FVFModel
(
    const word& name,
    const dictionary& phaseDict,
    const fvMesh& mesh
)
    :
    name_(name),
    phaseDict_(phaseDict),
    mesh_(mesh),
    rFVF_
    (
        IOobject
        (
            phaseDict.name()+".rFVF",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rFVF", dimless, 1.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    drFVFdP_
    (
        IOobject
        (
            phaseDict_.name()+".drFVFdP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("drFVFdP", dimless/dimPressure, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    )
{}

// ************************************************************************* //
