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
    Foam::samplePhase

Description

SourceFiles
    samplePhase.C

\*---------------------------------------------------------------------------*/

#ifndef samplePhase_H
#define samplePhase_H

#include "phase.H"
#include "fvcReconstruct.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phases
{

/*---------------------------------------------------------------------------*\
                         Class samplePhase Declaration
\*---------------------------------------------------------------------------*/

class samplePhase
:
    public phase
{
public:

    //- Runtime type information
    ClassName("samplePhase");


    // Constructors

        // Construct from components
        samplePhase
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const mixtureType& mT = mixtureType::singlePhase
        );

        // Construct from copy
        samplePhase (const samplePhase&);

	// Destructor
	virtual ~samplePhase();

    // Member Functions

		//- The actual interpolation operation
		virtual void correct()
        {
            U_ = fvc::reconstruct(phi());
            U_.correctBoundaryConditions();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interpolationTables

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "samplePhase.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
