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
    Foam::DiagAnisoRock

Description
    A template base class for aniostropic rocks where absoluet permeability
    is expressed as a normalized diagonal tensor.

SourceFiles
    DiagAnisoRock.C

\*---------------------------------------------------------------------------*/

#ifndef DiagAnisoRock_H
#define DiagAnisoRock_H

#include "UniformityTypes.H"
#include "IsotropyTypes.H"
#include "rock.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class DiagAnisoRock Declaration
\*---------------------------------------------------------------------------*/

template<class CompressiblityType>
class DiagAnisoRock
:
    public rock<DiagAnisotropic, CompressiblityType>
{
public:
    
    //- Runtime type information
    TypeName("DiagAnisoRock");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        DiagAnisoRock,
        dictionary,
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& rockProperties
        ),
        (name, mesh, rockProperties)
    );

    // Constructors

        // Construct from components
        DiagAnisoRock
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& rockProperties
        );

        // Construct from copy
        DiagAnisoRock (const DiagAnisoRock&);

    // Selectors

        //- Return a reference to the selected DiagAnisoRock
        static autoPtr<DiagAnisoRock> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& rockProperties
        );


    //- Destructor
    virtual ~DiagAnisoRock();

    // Member Operators

        //- Disallow default bitwise assignment
        DiagAnisoRock& operator=(const DiagAnisoRock&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DiagAnisoRock.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
