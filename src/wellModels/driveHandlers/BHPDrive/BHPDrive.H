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
    Foam::driveHandlers::BHPDrive

Description
    Imposes the Bottom-Hole-Pressure on a well by sending coeffs into the main
    matrices, and treating the pressure implicitely. Assuming BHP is constant
    in:

    \f$ q_i = a_i p_i + b_i BHP + c_i \f$

    This produces a diagonal matrix in p. 

SourceFiles
    BHPDrive.C

\*---------------------------------------------------------------------------*/

#ifndef BHPDrive_H
#define BHPDrive_H

#include "driveHandler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace driveHandlers
{

/*---------------------------------------------------------------------------*\
                        Class BHPDrive Declaration
\*---------------------------------------------------------------------------*/

template<class RockType, int nPhases>
class BHPDrive
:
    public driveHandler<RockType, nPhases>
{
protected:

    // Protected Data Members

        //- Attached phase
        wordList phases_;

public:

    //- Runtime type information
    ClassName("BHP");

    // Constructors

        // Construct from components
        BHPDrive
        (
            const word& name,
            const dictionary& driveDict,
            HashTable<autoPtr<wellSource<RockType, nPhases>>>& sources,
            sourceProperties& srcProps,
            HashPtrTable<fvScalarMatrix>& matrices
        );

        // Construct from copy
        BHPDrive(const BHPDrive&) = delete;

	// Destructor
	virtual ~BHPDrive();

    // Member Functions

        //- Update drive coeffs
        virtual void correct();

    // Member operators

        //- Disallow default bitwise assignment
        BHPDrive& operator=(const BHPDrive&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace driveHandlers

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BHPDrive.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
