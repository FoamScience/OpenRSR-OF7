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
    Foam::basicIncompressiblePhase

Description
    Standard Incompressible fluid classes possibly accounting for
    space-dependent viscosity


SourceFiles
    basicIncompressiblePhase.C

\*---------------------------------------------------------------------------*/

#ifndef basicIncompressiblePhase_H
#define basicIncompressiblePhase_H

#include "UniformityTypes.H"
#include "IncompressiblePhase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phases
{

/*---------------------------------------------------------------------------*\
                   Class basicIncompressiblePhase Declaration
\*---------------------------------------------------------------------------*/

template<class ViscosityType>
class basicIncompressiblePhase
:
    public IncompressiblePhase<ViscosityType>
{
public:

    //- Runtime type information
    ClassName("incompressible");


    // Constructors

        // Construct from components
        basicIncompressiblePhase
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const mixtureType& = mixtureType::singlePhase
        );

        // Construct from copy
        basicIncompressiblePhase (const basicIncompressiblePhase&);

	// Destructor
	virtual ~basicIncompressiblePhase();

    // Member Functions

		//- Phase fields update
		virtual void correct() override;

    // Member operators

        //- Disallow default bitwise assignment
        basicIncompressiblePhase& operator=(const basicIncompressiblePhase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phases

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicIncompressiblePhase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
