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

#include "VolatileDimensionedField.H"
#include "dimensionedType.H"
#include "Field.H"
#include "VolatileDimensionedFieldReuseFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
VolatileDimensionedField<Type, GeoMesh>::VolatileDimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& dims,
    const Field<Type>& field
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, dims, field),
    mask_(GeoMesh::size(mesh) == 1 ? 0 : -1)
{
}


template<class Type, class GeoMesh>
VolatileDimensionedField<Type, GeoMesh>::VolatileDimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& dims,
    const bool checkIOFlags
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, dims, checkIOFlags),
    mask_(GeoMesh::size(mesh) == 1 ? 0 : -1)
{
}


template<class Type, class GeoMesh>
VolatileDimensionedField<Type, GeoMesh>::VolatileDimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const bool checkIOFlags
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, dt, checkIOFlags),
    mask_(GeoMesh::size(mesh) == 1 ? 0 : -1)
{
}


template<class Type, class GeoMesh>
VolatileDimensionedField<Type, GeoMesh>::VolatileDimensionedField
(
    const VolatileDimensionedField<Type, GeoMesh>& df
)
:
    DimensionedField<Type, GeoMesh>(df),
    mask_(GeoMesh::size(df.mesh()) == 1 ? 0 : -1)
{}


template<class Type, class GeoMesh>
VolatileDimensionedField<Type, GeoMesh>::VolatileDimensionedField
(
    VolatileDimensionedField<Type, GeoMesh>& df,
    bool reuse
)
:
    DimensionedField<Type, GeoMesh>(df, reuse),
    mask_(GeoMesh::size(df.mesh()) == 1 ? 0 : -1)
{}


template<class Type, class GeoMesh>
VolatileDimensionedField<Type, GeoMesh>::VolatileDimensionedField
(
    VolatileDimensionedField<Type, GeoMesh>&& df
)
:
    DimensionedField<Type, GeoMesh>(df),
    mask_(move(df.mask_))
{}


template<class Type, class GeoMesh>
VolatileDimensionedField<Type, GeoMesh>::VolatileDimensionedField
(
    const tmp<VolatileDimensionedField<Type, GeoMesh>>& tdf
)
:
    DimensionedField<Type, GeoMesh>(tdf()),
    mask_(GeoMesh::size(this->mesh()) == 1 ? 0 : -1)
{
    tdf.clear();
}


template<class Type, class GeoMesh>
VolatileDimensionedField<Type, GeoMesh>::VolatileDimensionedField
(
    const IOobject& io,
    const VolatileDimensionedField<Type, GeoMesh>& df
)
:
    DimensionedField<Type, GeoMesh>(io, df),
    mask_(df.mask())
{}


template<class Type, class GeoMesh>
VolatileDimensionedField<Type, GeoMesh>::VolatileDimensionedField
(
    const IOobject& io,
    VolatileDimensionedField<Type, GeoMesh>& df,
    bool reuse
)
:
    DimensionedField<Type, GeoMesh>(io, df, reuse),
    mask_(df.mask())
{}


template<class Type, class GeoMesh>
VolatileDimensionedField<Type, GeoMesh>::VolatileDimensionedField
(
    const word& newName,
    const VolatileDimensionedField<Type, GeoMesh>& df
)
:
    DimensionedField<Type, GeoMesh>(newName, df),
    mask_(df.mask())
{}


template<class Type, class GeoMesh>
VolatileDimensionedField<Type, GeoMesh>::VolatileDimensionedField
(
    const word& newName,
    VolatileDimensionedField<Type, GeoMesh>& df,
    bool reuse
)
:
    DimensionedField<Type, GeoMesh>(newName, df, reuse),
    mask_(df.mask())
{}


template<class Type, class GeoMesh>
VolatileDimensionedField<Type, GeoMesh>::VolatileDimensionedField
(
    const word& newName,
    const tmp<VolatileDimensionedField<Type, GeoMesh>>& tdf
)
:
    DimensionedField<Type, GeoMesh>(newName, tdf()),
    mask_(GeoMesh::size(this->mesh()) == 1 ? 0 : -1)
{
    tdf.clear();
}


template<class Type, class GeoMesh>
tmp<VolatileDimensionedField<Type, GeoMesh>>
VolatileDimensionedField<Type, GeoMesh>::clone() const
{
    return tmp<VolatileDimensionedField<Type, GeoMesh>>
    (
        new VolatileDimensionedField<Type, GeoMesh>(*this)
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::VolatileDimensionedField<Type, GeoMesh>>
VolatileDimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& ds
)
{
    return tmp<VolatileDimensionedField<Type, GeoMesh>>
    (
        new VolatileDimensionedField<Type, GeoMesh>
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            ds,
            false
        )
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::VolatileDimensionedField<Type, GeoMesh>>
VolatileDimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensioned<Type>& dt
)
{
    return tmp<VolatileDimensionedField<Type, GeoMesh>>
    (
        new VolatileDimensionedField<Type, GeoMesh>
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dt,
            false
        )
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::VolatileDimensionedField<Type, GeoMesh>>
VolatileDimensionedField<Type, GeoMesh>::New
(
    const word& newName,
    const VolatileDimensionedField<Type, GeoMesh>& df
)
{
    return tmp<VolatileDimensionedField<Type, GeoMesh>>
    (
        new VolatileDimensionedField<Type, GeoMesh>
        (
            IOobject
            (
                newName,
                df.instance(),
                df.local(),
                df.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            df
        )
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::VolatileDimensionedField<Type, GeoMesh>>
VolatileDimensionedField<Type, GeoMesh>::New
(
    const word& newName,
    const tmp<VolatileDimensionedField<Type, GeoMesh>>& tdf
)
{
    return tmp<VolatileDimensionedField<Type, GeoMesh>>
    (
        new VolatileDimensionedField<Type, GeoMesh>
        (
            IOobject
            (
                newName,
                tdf().instance(),
                tdf().local(),
                tdf().db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            tdf
        )
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
VolatileDimensionedField<Type, GeoMesh>::~VolatileDimensionedField()
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

MEMBER_PRODUCT_OPERATOR(typeOfSum, +, Plus, add, VolatileDimensionedField)
MEMBER_PRODUCT_OPERATOR(typeOfSum, -, Minus,subtract, VolatileDimensionedField)
MEMBER_PRODUCT_OPERATOR(outerProduct, *, Star, outer, VolatileDimensionedField)
MEMBER_PRODUCT_OPERATOR(crossProduct, ^, Cross,cross, VolatileDimensionedField)
MEMBER_PRODUCT_OPERATOR(innerProduct, &, Inner, dot, VolatileDimensionedField)
MEMBER_PRODUCT_OPERATOR(scalarProduct,&&,Scalar,dotdot,VolatileDimensionedField)

TMP_PRODUCT_OPERATOR
(
    typeOfSum, +, Plus, add,
    VolatileDimensionedField, VolatileDimensionedField
)
TMP_PRODUCT_OPERATOR
(
    typeOfSum, -, Minus, subtract,
    VolatileDimensionedField, VolatileDimensionedField
)
TMP_PRODUCT_OPERATOR
(
    outerProduct, *, Star, outer,
    VolatileDimensionedField, VolatileDimensionedField
)
TMP_PRODUCT_OPERATOR
(
    crossProduct, ^, Cross, cross,
    VolatileDimensionedField, VolatileDimensionedField
)
TMP_PRODUCT_OPERATOR
(
    innerProduct, &, Inner, dot,
    VolatileDimensionedField, VolatileDimensionedField
)
TMP_PRODUCT_OPERATOR
(
    scalarProduct, &&, Scalar, dotdot,
    VolatileDimensionedField, VolatileDimensionedField
)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef checkField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VolatileDimensionedFieldIO.C"

// ************************************************************************* //
