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
    Foam::driveHandlers::flowRateDrive

Description
    Imposes the phase-flow-rate on a well by passing the requested value.

    \f$ q_i = a_i p_i + b_i BHP + c_i \f$

    Assuming \f$ \Sigma q_i = q_t\f$ is imposed,
    \f$ q_t = \Sigma (a_i p_i) + BHP \Sigma b_i + \Sigma c_i \f$
    gives a BHP value to be substituted in the cell's source equation

Note
    Treating pressure in neighbouring cells in the well implicitely.

SourceFiles
    flowRateDrive.C

\*---------------------------------------------------------------------------*/

#ifndef flowRateDrive_H
#define flowRateDrive_H

#include "driveHandler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace driveHandlers
{

/*---------------------------------------------------------------------------*\
                        Class flowRateDrive Declaration
\*---------------------------------------------------------------------------*/

template<class RockType, int nPhases>
class flowRateDrive
:
    public driveHandler<RockType, nPhases>
{
protected:

    // Protected Member Data

        //- Attached phase name
        word phase_;

public:

    //- Runtime type information
    ClassName("flowRate");

    // Constructors

        // Construct from components
        flowRateDrive
        (
            const word& name,
            const dictionary& driveDict,
            wellSource<RockType, nPhases>& source,
            sourceProperties& srcProps,
            HashPtrTable<fvScalarMatrix>& matrices
        );

        // Construct from copy
        flowRateDrive(const flowRateDrive&) = delete;

	// Destructor
	virtual ~flowRateDrive();

    // Member Functions

        //- Update drive coeffs
        virtual void correct();

    // Member operators

        //- Disallow default bitwise assignment
        flowRateDrive& operator=(const flowRateDrive&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace driveHandlers

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "flowRateDrive.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
