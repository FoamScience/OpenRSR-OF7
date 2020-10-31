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
    Foam::twoPhaseCFLMethods::CoatsNo

Description
    Two-Phase Coats number calculation for porous mediums

Note
    For use in IMPES simulations

SourceFiles
    CoatsNo.C

\*---------------------------------------------------------------------------*/

#ifndef CoatsNo_H
#define CoatsNo_H

#include "phase.H"
#include "CFLMethod.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace twoPhaseCFLMethods
{

/*---------------------------------------------------------------------------*\
                        Class CoatsNo Declaration
\*---------------------------------------------------------------------------*/

template<class RockType>
class CoatsNo
:
    public CFLMethod<RockType, 2>
{
protected:

    // Protected Data Members

        //- Access to canonical and noncanonical phases
        wordList phaseNames_;
        const phase& cPhase_;
        const phase& nPhase_;

        //- Dictionary of default needed field names
        dictionary fieldNames_;

        //- Viscosity ratio
        volScalarField muRatio_;

        //- Flux differential
        scalarField dPhi_;

public:

    //- Runtime type information
    ClassName("CoatsNo");

    // Constructors

        // Construct from components
        CoatsNo
        (
            const word& name,
            const dictionary& algorithmProperties,
            const surfaceScalarField& phi,
            const wordList& phaseNames,
            const RockType& rock
        );

        // Construct from copy
        CoatsNo(const CoatsNo&) = delete;

	// Destructor
	virtual ~CoatsNo();

    // Member Functions

		//- Correct CFL number
		virtual void correct();

    // Member operators

        //- Disallow default bitwise assignment
        CoatsNo& operator=(const CoatsNo&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseCFLMethods

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CoatsNo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
