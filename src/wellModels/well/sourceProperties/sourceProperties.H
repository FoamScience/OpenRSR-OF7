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
    Foam::sourceProperties

Description
    A class to hold source-related well properties.

SourceFiles
    sourceProperties.C

\*---------------------------------------------------------------------------*/

#ifndef sourceProperties_H
#define sourceProperties_H

#include "fvMesh.H"
#include "dimensionedScalar.H"
#include "UniformDimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sourceProperties Declaration
\*---------------------------------------------------------------------------*/

class sourceProperties
{

public:

        //- Enum for handling well orientation
        enum class orientationHandling
        {
            // TODO: Implement other directions
            vertical, generic
        };

protected:

    // Protected Data Members
    
        //- Gravitational acceleration
        const UniformDimensionedField<vector>& g_;

        //- well orientation
        orientationHandling orientation_;

        //- Well radius
        dimensionedScalar radius_;

        //- Well skin
        scalar skin_;

        //- Well index
        scalarList J_;

public:

    //- Runtime type information
    ClassName("sourceProperties");

    // Constructors 

        //- Construct from dictionary and mesh
        sourceProperties(const fvMesh& mesh, const dictionary& wellDict_);

    //- Destructor
    virtual ~sourceProperties() {}

    //- Return orientation handling as word
    static word orientationHandlingToWord(const orientationHandling& ori);

    //- Return orientation handling as enumeration
    static orientationHandling wordToOrientationHandling(const word& ori);

    //- Return orientation
    const orientationHandling& orientation() const
    {
        return orientation_;
    }

    //- Return gravitational acceleration
    const UniformDimensionedField<vector>& g() const
    {
        return g_;
    }

    //- Return radius
    const dimensionedScalar& radius() const
    {
        return radius_;
    }

    //- Return skin
    const scalar& skin() const
    {
        return skin_;
    }

    //- Return well index
    const scalarList& wellIndex() const
    {
        return J_;
    }

    //- Access to well index
    scalarList& wellIndex() {
        return J_;
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
