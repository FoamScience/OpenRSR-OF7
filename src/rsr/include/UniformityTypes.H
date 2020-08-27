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
    Uniformity Types for scalar values
    
    - Compressiblity: Compressible, Incompressible

\*---------------------------------------------------------------------------*/

#ifndef UniformityTypes_H
#define UniformityTypes_H

#include "UniformDimensionedField.H"
#include "volFields.H"
#include "fvMesh.H"
#include "calculatedFvPatchFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Compressible Scalar Type
using Compressible = volScalarField;

//- Incompressible Scalar Type
template<class cmptType> class UniformVolumeField;
using Incompressible = UniformVolumeField<scalar>;

//- Incompressible Type
template<class cmptType>
class UniformVolumeField
:
    public UniformDimensionedField<cmptType>
{

public :

    // Constructors

        //- Construct from Compressible Components
        UniformVolumeField
        (
            const IOobject& io,
            const fvMesh& mesh,
            const dimensioned<cmptType>& dt,
            const word& patchFieldType = fvPatchField<cmptType>::calculatedType()
        )
        : UniformDimensionedField<cmptType>(io, dt)
        {}

        //- Construct from copy
        UniformVolumeField(const UniformVolumeField& incomp)
        :
            UniformDimensionedField<cmptType>(incomp)
        {
        }

        //- Construct from Istream
        UniformVolumeField(const IOobject& io)
        :
            UniformDimensionedField<cmptType>(io)
        {
        }

        //- Move constructor
        UniformVolumeField(const UniformVolumeField&&) = delete;

    //- Destructor
    virtual ~UniformVolumeField(){}

    //- Member operators

        //- Disallow default bitwise assignment
        UniformVolumeField<cmptType>& operator=(const UniformVolumeField& other)
        {
            *this = UniformDimensionedField<cmptType>::operator=(other);
            return *this;
        }

        //- Move operator
        UniformVolumeField<cmptType>&
        operator=(UniformVolumeField&&) noexcept = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
