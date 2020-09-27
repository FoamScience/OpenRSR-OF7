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

InClass
    Foam::VolatileDimensionedField

Description
    Macros for declaring and defining operators and functions for volatile
    field operations

\*---------------------------------------------------------------------------*/

#ifndef VolatileDimensionedFieldM_H
#define VolatileDimensionedFieldM_H

#include "products.H"
#include "VolatileDimensionedFieldOpsTraitsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Check mesh for two volatile fields
// require at least a  single-cell meshe if meshes are different
#define checkField(df1, df2, op)                                               \
if                                                                             \
(                                                                              \
    ( GeoMesh::size((df1).mesh()) != 1 )                                       \
    and ( GeoMesh::size((df2).mesh()) != 1 )                                   \
    and ( &(df1).mesh() != &(df2).mesh() )                                     \
)                                                                              \
{                                                                              \
    FatalErrorInFunction                                                       \
        << "Different mesh for fields "                                        \
        << (df1).name() << " and " << (df2).name()                             \
        << " during operation "  #op                                           \
        << abort(FatalError);                                                  \
}

#define MEMBER_PRODUCT_OPERATOR_DECL(product, op, opFunc, fieldType)           \
                                                                               \
template<class Type2>                                                          \
tmp<VolatileDimensionedField<typename product<Type, Type2>::type, GeoMesh>>    \
operator op                                                                    \
(                                                                              \
    const fieldType<Type2, GeoMesh>& df2                                       \
);                                                                             \
                                                                               \
template<class Type2>                                                          \
tmp<VolatileDimensionedField<typename product<Type, Type2>::type, GeoMesh>>    \
operator op                                                                    \
(                                                                              \
    const tmp<fieldType<Type2, GeoMesh>>& tdf2                                 \
)

#define MEMBER_PRODUCT_OPERATOR(product, op, opName, opFunc, fieldType)        \
                                                                               \
BEST_BINARY_PRODUCT(op, opName, opFunc, );                                     \
                                                                               \
template<class Type, class GeoMesh>                                            \
template<class Type2>                                                          \
tmp<VolatileDimensionedField<typename product<Type, Type2>::type, GeoMesh>>    \
VolatileDimensionedField<Type, GeoMesh>::operator op                           \
(                                                                              \
    const fieldType<Type2, GeoMesh>& df2                                       \
)                                                                              \
{                                                                              \
    checkField((*this), df2, op);                                              \
    typedef typename product<Type, Type2>::type productType;                   \
    tmp<VolatileDimensionedField<productType, GeoMesh>> tRes                   \
    (                                                                          \
        VolatileDimensionedField<productType, GeoMesh>::New                    \
        (                                                                      \
            '(' + this->name() + #op + df2.name() + ')',                       \
            this->mask() ? this->mesh() : df2.mesh(),                          \
            this->dimensions() op df2.dimensions()                             \
        )                                                                      \
    );                                                                         \
                                                                               \
    forAll(tRes(), ci)                                                         \
    {                                                                          \
        opName##Helper::opFunc##Best                                           \
        (                                                                      \
            tRes.ref()[ci], this->operator[](ci), df2[ci]                      \
        );                                                                     \
    }                                                                          \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type, class GeoMesh>                                            \
template<class Type2>                                                          \
tmp<VolatileDimensionedField<typename product<Type, Type2>::type, GeoMesh>>    \
VolatileDimensionedField<Type, GeoMesh>::operator op                           \
(                                                                              \
    const tmp<fieldType<Type2, GeoMesh>>& tdf2                                 \
)                                                                              \
{                                                                              \
    checkField((*this), tdf2(), op);                                           \
    typedef typename product<Type, Type2>::type productType;                   \
                                                                               \
    const VolatileDimensionedField<Type2, GeoMesh>& df2 = tdf2();              \
                                                                               \
    tmp<VolatileDimensionedField<productType, GeoMesh>> tRes =                 \
        reuseTmpVolatileDimensionedField<productType, Type2, GeoMesh>::New     \
        (                                                                      \
            tdf2,                                                              \
            '(' + this->name() + #op + df2.name() + ')',                       \
            this->dimensions() op df2.dimensions()                             \
        );                                                                     \
                                                                               \
    forAll(tRes(), ci)                                                         \
    {                                                                          \
        opName##Helper::opFunc##Best                                           \
        (                                                                      \
            tRes.ref()[ci], this->operator[](ci), df2[ci]                      \
        );                                                                     \
    }                                                                          \
                                                                               \
    tdf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                               

#define TMP_PRODUCT_OPERATOR_DECL\
(product, op, opName, opFunc, fieldType1, fieldType2)                          \
                                                                               \
template<class Type1, class Type2, class GeoMesh>                              \
tmp<VolatileDimensionedField<typename product<Type1, Type2>::type, GeoMesh>>   \
operator op                                                                    \
(                                                                              \
    const tmp<fieldType1<Type1, GeoMesh>>& tdf1,                               \
    const tmp<fieldType2<Type2, GeoMesh>>& tdf2                                \
)                                                                           

#define TMP_PRODUCT_OPERATOR\
(product, op, opName, opFunc, fieldType1, fieldType2)                          \
                                                                               \
BEST_BINARY_PRODUCT(op, opName, opFunc, Tmp);                                  \
                                                                               \
template<class Type1, class Type2, class GeoMesh>                              \
tmp<VolatileDimensionedField<typename product<Type1, Type2>::type, GeoMesh>>   \
operator op                                                                    \
(                                                                              \
    const tmp<fieldType1<Type1, GeoMesh>>& tdf1,                               \
    const tmp<fieldType2<Type2, GeoMesh>>& tdf2                                \
)                                                                              \
{                                                                              \
    checkField(tdf1(), tdf2(), op);                                            \
    typedef typename product<Type1, Type2>::type productType;                  \
                                                                               \
    const VolatileDimensionedField<Type1, GeoMesh>& df1 = tdf1();              \
    const VolatileDimensionedField<Type2, GeoMesh>& df2 = tdf2();              \
                                                                               \
    tmp<VolatileDimensionedField<productType, GeoMesh>> tRes =                 \
        reuseTmpTmpVolatileDimensionedField                                    \
        <productType, Type1, Type1, Type2, GeoMesh>::New                       \
        (                                                                      \
            !df1.mask() ? tdf1 : tdf2,                                         \
            df1.mask() ? tdf2 : tdf1,                                          \
            '(' + df1.name() + #op + df2.name() + ')',                         \
            df1.dimensions() op df2.dimensions()                               \
        );                                                                     \
                                                                               \
    forAll(tRes(), ci)                                                         \
    {                                                                          \
        opName##HelperTmp::opFunc##Best                                        \
        (                                                                      \
            tRes.ref()[ci], df1[ci], df2[ci]                                   \
        );                                                                     \
    }                                                                          \
                                                                               \
    tdf1.clear();                                                              \
    tdf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
