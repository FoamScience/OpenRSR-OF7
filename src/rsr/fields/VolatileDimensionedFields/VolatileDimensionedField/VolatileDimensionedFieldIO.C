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

#include "IOstreams.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::VolatileDimensionedField<Type, GeoMesh>::VolatileDimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const word& fieldDictEntry
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, fieldDictEntry),
    mask_(GeoMesh::size(mesh) == 1 ? 0 : -1)
{
}


template<class Type, class GeoMesh>
Foam::VolatileDimensionedField<Type, GeoMesh>::VolatileDimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dictionary& fieldDict,
    const word& fieldDictEntry
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, fieldDictEntry),
    mask_(GeoMesh::size(mesh) == 1 ? 0 : -1)
{
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const VolatileDimensionedField<Type, GeoMesh>& df
)
{
    df.writeData(os);

    return os;
}


template<class Type, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tmp<VolatileDimensionedField<Type, GeoMesh>>& tdf
)
{
    tdf().writeData(os);
    tdf.clear();

    return os;
}


// ************************************************************************* //
