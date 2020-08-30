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
    Foam::basicIsotropicRock

Description
    Standard implementation of an isotropic rock class

SourceFiles
    basicIsotropicRock.C

\*---------------------------------------------------------------------------*/

#ifndef basicIsotropicRock_H
#define basicIsotropicRock_H

#include "IsoRock.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace rocks
{

/*---------------------------------------------------------------------------*\
                   Class basicIsotropicRock Declaration
\*---------------------------------------------------------------------------*/

template<class CompressibilityType>
class basicIsotropicRock
:
    public IsoRock<CompressibilityType>
{
public:

    //- Runtime type information
    ClassName("isotropic");


    // Constructors

        // Construct from components
        basicIsotropicRock
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& rockProperties
        );

        // Construct from copy
        basicIsotropicRock (const basicIsotropicRock&);

	// Destructor
	virtual ~basicIsotropicRock();

    // Member Functions

		//- Phase fields update
		virtual void correct() override;

    // Member operators

        //- Disallow default bitwise assignment
        basicIsotropicRock& operator=(const basicIsotropicRock&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phases

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicIsotropicRock.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
