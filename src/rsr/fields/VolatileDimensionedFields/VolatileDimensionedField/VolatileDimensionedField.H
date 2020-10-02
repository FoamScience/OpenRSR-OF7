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

Class
    Foam::VolatileDimensionedField

Description
    A field with dimensions. Can be instantiated on a single-cell mesh or
    on a regular volume mesh. Mainly for optimized operations between
    dimensioned fields and uniform-like ones.

SourceFiles
    VolatileDimensionedFieldI.H
    VolatileDimensionedField.C
    VolatileDimensionedFieldIO.C

\*---------------------------------------------------------------------------*/

#ifndef VolatileDimensionedField_H
#define VolatileDimensionedField_H

#include "DimensionedField.H"
#include "surfaceMesh.H"
#include "VolatileDimensionedFieldM.H"
#include "products.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type, class GeoMesh> class VolatileDimensionedField;

template<class Type, class GeoMesh> Ostream& operator<<
(
    Ostream&,
    const VolatileDimensionedField<Type, GeoMesh>&
);

template<class Type, class GeoMesh> Ostream& operator<<
(
    Ostream&,
    const tmp<VolatileDimensionedField<Type, GeoMesh>>&
);

/*---------------------------------------------------------------------------*\
                           Class VolatileDimensionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class VolatileDimensionedField
:
    public DimensionedField<Type, GeoMesh>
{

public:

    // Public Typedefs

        //- Type of mesh on which this VolatileDimensionedField is instantiated
        typedef typename GeoMesh::Mesh Mesh;

        //- Type of the field from which this VolatileDimensionedField is derived
        typedef Field<Type> FieldType;

        //- Component type of the elements of the field
        typedef typename Field<Type>::cmptType cmptType;


private:

    // Private Data

        //- A mask to decide on how much space to allocate for field values
        label mask_;

public:

    //- Runtime type information
    ClassName("VolatileDimensionedField");

    // Constructors

        //- Construct from components
        VolatileDimensionedField
        (
            const IOobject&,
            const Mesh& mesh,
            const dimensionSet&,
            const Field<Type>&
        );

        //- Construct from components
        //  Used for temporary fields which are initialised after construction
        VolatileDimensionedField
        (
            const IOobject&,
            const Mesh& mesh,
            const dimensionSet&,
            const bool checkIOFlags = true
        );

        //- Construct from components
        VolatileDimensionedField
        (
            const IOobject&,
            const Mesh& mesh,
            const dimensioned<Type>&,
            const bool checkIOFlags = true
        );

        //- Construct from Istream
        VolatileDimensionedField
        (
            const IOobject&,
            const Mesh& mesh,
            const word& fieldDictEntry="value"
        );

        //- Construct from dictionary
        VolatileDimensionedField
        (
            const IOobject&,
            const Mesh& mesh,
            const dictionary& fieldDict,
            const word& fieldDictEntry="value"
        );

        //- Copy constructor
        VolatileDimensionedField
        (
            const VolatileDimensionedField<Type, GeoMesh>&
        );

        //- Copy constructor or re-use as specified.
        VolatileDimensionedField
        (
            VolatileDimensionedField<Type, GeoMesh>&,
            bool reuse
        );

        //- Move constructor
        VolatileDimensionedField
        (
            VolatileDimensionedField<Type, GeoMesh>&&
        );

        //- Copy constructor of tmp<VolatileDimensionedField> deleting argument
        VolatileDimensionedField
        (
            const tmp<VolatileDimensionedField<Type, GeoMesh>>&
        );

        //- Copy constructor resetting IO parameters
        VolatileDimensionedField
        (
            const IOobject&,
            const VolatileDimensionedField<Type, GeoMesh>&
        );

        //- Copy constructor resetting IO parameters and re-use as specified.
        VolatileDimensionedField
        (
            const IOobject&,
            VolatileDimensionedField<Type, GeoMesh>&,
            bool reuse
        );

        //- Copy constructor resetting name
        VolatileDimensionedField
        (
            const word& newName,
            const VolatileDimensionedField<Type, GeoMesh>&
        );

        //- Copy constructor resetting name and re-use as specified.
        VolatileDimensionedField
        (
            const word& newName,
            VolatileDimensionedField<Type, GeoMesh>&,
            bool reuse
        );

        //- Copy constructor resetting name
        VolatileDimensionedField
        (
            const word& newName,
            const tmp<VolatileDimensionedField<Type, GeoMesh>>&
        );

        //- Clone
        tmp<VolatileDimensionedField<Type, GeoMesh>> clone() const;

        //- Return a temporary field constructed from name, mesh
        //  and dimensionSet
        static tmp<VolatileDimensionedField<Type, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensionSet&
        );

        //- Return a temporary field constructed from name, mesh
        //  and dimensionType providing dimensions and value
        static tmp<VolatileDimensionedField<Type, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensioned<Type>&
        );

        //- Return a temporary field constructed from name and a field
        static tmp<VolatileDimensionedField<Type, GeoMesh>> New
        (
            const word& newName,
            const VolatileDimensionedField<Type, GeoMesh>&
        );

        //- Return renamed temporary field
        static tmp<VolatileDimensionedField<Type, GeoMesh>> New
        (
            const word& newName,
            const tmp<VolatileDimensionedField<Type, GeoMesh>>&
        );


    //- Destructor
    virtual ~VolatileDimensionedField();


    // Member Functions

        //- Return mask
        inline const label& mask() const
        {
            return mask_;
        }

    // Member Operators

        // Access to field elements
        Type& operator[](const label& index)
        {
            return DimensionedField<Type, GeoMesh>::operator[]( index & mask_);
        }

        // Const-access to field elements
        const Type& operator[](const label& index) const
        {
            return DimensionedField<Type, GeoMesh>::operator[]( index & mask_);
        }

        // Product operators using operator between field elements
        MEMBER_PRODUCT_OPERATOR_DECL
        (
            typeOfSum, +, add, VolatileDimensionedField
        );
        MEMBER_PRODUCT_OPERATOR_DECL
        (
            typeOfSum, -, subtract, VolatileDimensionedField
        );
        MEMBER_PRODUCT_OPERATOR_DECL
        (
            outerProduct, *, outer, VolatileDimensionedField
        );
        MEMBER_PRODUCT_OPERATOR_DECL
        (
            crossProduct, ^, cross, VolatileDimensionedField
        );
        MEMBER_PRODUCT_OPERATOR_DECL
        (
            innerProduct, &, dot, VolatileDimensionedField
        );
        MEMBER_PRODUCT_OPERATOR_DECL
        (
            scalarProduct, &&, dotdot, VolatileDimensionedField
        );


    // Ostream Operators

        friend Ostream& operator<< <Type, GeoMesh>
        (
            Ostream&,
            const VolatileDimensionedField<Type, GeoMesh>&
        );

        friend Ostream& operator<< <Type, GeoMesh>
        (
            Ostream&,
            const tmp<VolatileDimensionedField<Type, GeoMesh>>&
        );

};

TMP_PRODUCT_OPERATOR_DECL
(
    typeOfSum, +, Plus, add,
    VolatileDimensionedField, VolatileDimensionedField
);
TMP_PRODUCT_OPERATOR_DECL
(
    typeOfSum, -, Minus, subtract,
    VolatileDimensionedField, VolatileDimensionedField
);
TMP_PRODUCT_OPERATOR_DECL
(
    outerProduct, *, Star, outer,
    VolatileDimensionedField, VolatileDimensionedField
);
TMP_PRODUCT_OPERATOR_DECL
(
    crossProduct, ^, Cross, cross,
    VolatileDimensionedField, VolatileDimensionedField
);
TMP_PRODUCT_OPERATOR_DECL
(
    innerProduct, &, Inner, dot,
    VolatileDimensionedField, VolatileDimensionedField
);
TMP_PRODUCT_OPERATOR_DECL
(
    scalarProduct, &&, Scalar, dotdot,
    VolatileDimensionedField, VolatileDimensionedField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "VolatileDimensionedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
