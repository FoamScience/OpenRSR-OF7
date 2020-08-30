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

#include "DiagAnisoRock.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CompressibilityType>
Foam::autoPtr<Foam::DiagAnisoRock<CompressibilityType>>
Foam::DiagAnisoRock<CompressibilityType>::New
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& rockProperties
)
{
    const word rockType = rockProperties.subDict(name).lookupOrDefault<word>
    (
      "rockType",
      "diagAnisotropic"
    );

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(rockType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown DiagAnisotropic Rock type "
            << rockType << nl << nl
            << "Valid DiagAnisotropic Rock types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<DiagAnisoRock>
    ( cstrIter()(name, mesh, rockProperties) );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompressibilityType>
Foam::DiagAnisoRock<CompressibilityType>::DiagAnisoRock
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& rockProperties
)
:
    rock<DiagAnisotropic, CompressibilityType>(name, mesh, rockProperties)
{
}


template<class CompressibilityType>
Foam::DiagAnisoRock<CompressibilityType>::DiagAnisoRock
(
    const DiagAnisoRock& rk
)
:
    rock<DiagAnisotropic, CompressibilityType>(rk)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompressibilityType>
Foam::DiagAnisoRock<CompressibilityType>::~DiagAnisoRock()
{}

// ************************************************************************* //
