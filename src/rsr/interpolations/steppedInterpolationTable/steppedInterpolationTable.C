/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "steppedInterpolationTable.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationTables::steppedInterpolationTable<Type>::
steppedInterpolationTable
(
    const List<Tuple2<scalar, List<Type>>>& values,
	const bool isPeriodic
)
:
    basicInterpolationTable<Type>(values, isPeriodic)
{
}


template<class Type>
Foam::interpolationTables::steppedInterpolationTable<Type>::
steppedInterpolationTable
(
    const dictionary& dict
)
:
    basicInterpolationTable<Type>(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationTables::steppedInterpolationTable<Type>::
~steppedInterpolationTable()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::List<Type> 
Foam::interpolationTables::steppedInterpolationTable<Type>::interpolate
(
    const Foam::scalar& time
) const
{
    // Always return the single value if data size == 1
	if (this->values_.size() <= 1)
	{
		return this->values_[0].second();
	}
	
    // Correct scalar value if periodicity is enabled
    Foam::scalar pTime = this->projectTime(time);

    // Find index to next (or lookup) element in values list
	int nextElement = this->lookup(pTime);

    // Extract upstream element from values list
	const Tuple2<scalar, List<Type>>& us = this->values_[nextElement];
	if(us.first() == pTime)
	{
		return us.second();
	}

    // nextElement can't be 0
	return this->values_[nextElement-1].second();
}

// ************************************************************************* //
