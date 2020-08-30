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
    IsotropicDimentionedField

Description
    Isotropy/Anisotropy Types with identical construction interface
    
    - Isotropic: Identical field in all directions
    - Anisotropic: Tensorial representation of field changes per direction
    - DiagAnisotropic: Diagonal tensor representation of field changes
                       per main-axis direction

\*---------------------------------------------------------------------------*/

#ifndef IsotropyTypes_H
#define IsotropyTypes_H

#include "UniformDimensionedField.H"
#include "pTraits.H"
#include "volFields.H"
#include "fvMesh.H"
#include "uniformDimensionedFields.H"
#include "volFieldsFwd.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Isotropic
using Isotropic = volScalarField;

//- DiagAnisotropic
using DiagAnisotropic = volVectorField;

//- Anisotropic
using Anisotropic = volTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
