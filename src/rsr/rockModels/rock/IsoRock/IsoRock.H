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
    Foam::IsoRock

Description
    A template base class for iostropic rocks.

SourceFiles
    IsoRock.C

\*---------------------------------------------------------------------------*/

#ifndef IsoRock_H
#define IsoRock_H

#include "UniformityTypes.H"
#include "IsotropyTypes.H"
#include "rock.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IsoRock Declaration
\*---------------------------------------------------------------------------*/

template<class CompressiblityType>
class IsoRock
:
    public rock<Isotropic, CompressiblityType>
{
public:
    
    //- Runtime type information
    TypeName("IsoRock");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        IsoRock,
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
        IsoRock
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& rockProperties
        );

        // Construct from copy
        IsoRock (const IsoRock&);

    // Selectors

        //- Return a reference to the selected IsoRock
        static autoPtr<IsoRock> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& rockProperties
        );


    //- Destructor
    virtual ~IsoRock();

    // Member Operators

        //- Disallow default bitwise assignment
        IsoRock& operator=(const IsoRock&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IsoRock.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
