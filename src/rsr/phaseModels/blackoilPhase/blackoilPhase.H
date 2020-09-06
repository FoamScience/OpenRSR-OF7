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
    Foam::blackoilPhase

Description
    Compressible fluid class, where compressibility is expressed in terms
    of FVF. Accounts for variant viscosity.


SourceFiles
    blackoilPhase.C

\*---------------------------------------------------------------------------*/

#ifndef blackoilPhase_H
#define blackoilPhase_H

#include "phase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phases
{

/*---------------------------------------------------------------------------*\
                        Class blackoilPhase Declaration
\*---------------------------------------------------------------------------*/

template<class ViscosityType>
class blackoilPhase
:
    public phase<Compressible, ViscosityType>
{
public:

    //- Runtime type information
    ClassName("blackoil");


    // Constructors

        // Construct from components
        blackoilPhase
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const mixtureType& = mixtureType::singlePhase
        );

        // Construct from copy
        blackoilPhase(const blackoilPhase&);

	// Destructor
	virtual ~blackoilPhase();

    // Member Functions

		//- Phase fields update
		virtual void correct() override;

    // Member operators

        //- Disallow default bitwise assignment
        blackoilPhase& operator=(const blackoilPhase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phases

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "blackoilPhase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
