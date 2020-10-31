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
    Foam::twoPhaseRelPermModels::krTabular

Description
    Two-Phase relative permeability model reading kr values and its derivatives
    as a function of canonical saturation from a table.

SourceFiles
    krTabular.C

\*---------------------------------------------------------------------------*/

#ifndef krTabular_H
#define krTabular_H

#include "dimensionedScalarFwd.H"
#include "basicInterpolationTable.H"
#include "relPermModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace twoPhaseRelPermModels
{

/*---------------------------------------------------------------------------*\
                        Class krTabular Declaration
\*---------------------------------------------------------------------------*/

template<class RockType>
class krTabular
:
    public relPermModel<RockType, 2>
{
    //- Name of the other phase
    word otherPhase_;

    //- Const-access to canonical saturation
    const volScalarField& alpha_;

    //- Time series of (Kr and dKr) vs canonical saturation
    autoPtr<basicInterpolationTable<scalar>> krSeries_;

public:

    //- Runtime type information
    ClassName("tabular");

    // Constructors

        // Construct from components
        krTabular
        (
            const word& name,
            const dictionary& transportProperties,
            const RockType& rock
        );

        // Construct from copy
        krTabular(const krTabular&) = delete;

	// Destructor
	virtual ~krTabular();

    // Member Functions

		//- Correct relative-permeability fields
		virtual void correct();

    // Member operators

        //- Disallow default bitwise assignment
        krTabular& operator=(const krTabular&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseRelPermModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "krTabular.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
