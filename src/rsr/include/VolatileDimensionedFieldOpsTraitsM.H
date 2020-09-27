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

Description
    Macros to pick the correct operations on VolatileDimensionedField
    components.

    If the expression (res = x op y) compiles, it's the preferred path.
    Else, we use opFunc(res, x, y) if that compiles.
    If nothing works, then, we don't know what to do and the compiler should
    ERROR OUT.

    Eg. The simple (res = x op y) may not compile because of unambiguous call
    to 'operator op' so we fall back to the less common opFunc(res, x, y)

Note
    Particularly, we don't want to fall back to any other form of the above
    operation because we can't trust it to "binary". Hence the compilation
    error instead of looking for other alternatives.

\*---------------------------------------------------------------------------*/

#ifndef VolatileDimensionedFieldOpsTraitsM_H
#define VolatileDimensionedFieldOpsTraitsM_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BEST_BINARY_PRODUCT(op, opName, opFunc, helperSuffix)                  \
                                                                               \
struct opName##Helper##helperSuffix                                            \
{                                                                              \
    /* See if (result = x op y) compiles */                                    \
    template<class rT, class T1, class T2>                                     \
    static auto operator##opName##Works(rT& res, T1& x ,T2& y)                 \
            -> decltype(res = x op y, std::true_type());                       \
    static std::false_type operator##opName##Works(...);                       \
                                                                               \
    /* See if add(result, x, y) compiles */                                    \
    template<class rT, class T1, class T2>                                     \
    static auto opFunc##Works(rT& res, T1& x , T2& y)                          \
    -> decltype(add(res, x, y), std::true_type());                             \
    static std::false_type opFunc##Works(...);                                 \
                                                                               \
    /* Selector for best path to take */                                       \
    template<class rT, class T1, class T2, bool has_op, bool has_func>         \
    struct best##opName##Selector;                                             \
                                                                               \
    /* Use this if (res = x op y) compiles */                                  \
    template<class rT, class T1, class T2, bool has_func>                      \
    struct best##opName##Selector<rT, T1, T2, true , has_func>                 \
    {                                                                          \
        static auto f(rT& res, T1&x, T2& y) ->decltype(res = x op y)           \
        {                                                                      \
            /*Info << "x " #op " y is happening" << endl;*/                    \
            res = x op y;                                                      \
            return res;                                                        \
        }                                                                      \
    };                                                                         \
                                                                               \
    /* Use this if opFunc(res,x,y) compiles */                                 \
    template<class rT, class T1, class T2>                                     \
    struct best##opName##Selector<rT, T1, T2, false , true>                    \
    {                                                                          \
        static auto f(rT& res, T1&x, T2& y) -> decltype(opFunc(res, x, y))     \
        {                                                                      \
            /*Info << "opFunc(res, x, y) is happening" << endl;*/              \
            return opFunc(res, x, y);                                          \
        }                                                                      \
    };                                                                         \
                                                                               \
    /* Compiler error if nothing works */                                      \
    template<class rT, class T1, class T2>                                     \
    struct best##opName##Selector<rT, T1, T2, false , false>                   \
    {                                                                          \
        static auto f(rT& res, T1&x, T2& y)                                    \
        {                                                                      \
            static_assert                                                      \
            (                                                                  \
                sizeof(rT) < 0,                                                \
                "No overload is appropriate for operation " #opName            \
            );                                                                 \
            return 0;                                                          \
        }                                                                      \
    };                                                                         \
                                                                               \
    /* Call this function to automatically determine what to do */             \
    template<class rT, class T1, class T2>                                     \
    static auto opFunc##Best(rT& res, T1& x, T2& y)                            \
    {                                                                          \
        return best##opName##Selector                                          \
            <                                                                  \
                rT, T1, T2,                                                    \
                decltype(operator##opName##Works(res, x, y))::value,           \
                decltype(opFunc##Works(res, x, y))::value                      \
            >::f(res, x, y);                                                   \
    }                                                                          \
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
