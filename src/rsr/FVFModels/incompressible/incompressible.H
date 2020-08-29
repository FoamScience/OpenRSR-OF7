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
    Foam::FVFModels::incompressible

Description
    Tabular description of the relationship 1|FVF = f(P). Typically used
    for undersaturated isotherm reservoirs.

SourceFiles
    incompressible.C

\*---------------------------------------------------------------------------*/

#ifndef incompressible_H
#define incompressible_H

#include "FVFModel.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FVFModels
{

/*---------------------------------------------------------------------------*\
                     Class incompressible Declaration
\*---------------------------------------------------------------------------*/

class incompressible
:
    public FVFModel<Incompressible>
{

public:

    //- Runtime type information
    ClassName("incompressible");

    // Constructors

        //- Construct from components
        incompressible
        (
            const word& name,
            const dictionary& phaseDict,
            const fvMesh& mesh
        );

        //- Construct from copy
        incompressible(const incompressible&);

    // Destructor
    virtual ~incompressible() {}

    // Member Functions

        //- Correct the FVF
        virtual void correct() override;

    // Member Operators

        //- Disallow default bitwise assignment
        incompressible& operator=(const incompressible&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FVFModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
