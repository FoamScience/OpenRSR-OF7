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
    Foam::basicDiagAnisotropicRock

Description
    Standard implementation of anisotropic rock class where permeability is
    expressed as a diagonal tensor

SourceFiles
    basicDiagAnisotropicRock.C

\*---------------------------------------------------------------------------*/

#ifndef basicDiagAnisotropicRock_H
#define basicDiagAnisotropicRock_H

#include "rock.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace rocks
{

/*---------------------------------------------------------------------------*\
                   Class basicDiagAnisotropicRock Declaration
\*---------------------------------------------------------------------------*/

class basicDiagAnisotropicRock
:
    public rock<DiagAnisotropic>
{
public:

    //- Runtime type information
    ClassName("standard");


    // Constructors

        // Construct from components
        basicDiagAnisotropicRock
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& rockProperties
        );

        // Construct from copy
        basicDiagAnisotropicRock (const basicDiagAnisotropicRock&);

	// Destructor
	virtual ~basicDiagAnisotropicRock();

    // Member Functions

		//- Phase fields update
		virtual void correct() override;

    // Member operators

        //- Disallow default bitwise assignment
        basicDiagAnisotropicRock& operator=
        (
            const basicDiagAnisotropicRock&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phases

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "basicDiagAnisotropicRock.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
