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

#include "interpolations/basicInterpolationTable/basicInterpolationTable.H"
#include "linearInterpolationTable.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationTables::linearInterpolationTable<Type>::
linearInterpolationTable
(
    const List<Tuple2<scalar, List<Type>>>& values,
	const bool isPeriodic
)
:
    basicInterpolationTable<Type>(values, isPeriodic)
{
}


template<class Type>
Foam::interpolationTables::linearInterpolationTable<Type>::
linearInterpolationTable
(
    const dictionary& dict
)
:
    basicInterpolationTable<Type>(dict)
{}


template<class Type>
Foam::interpolationTables::linearInterpolationTable<Type>::
linearInterpolationTable
(
    const linearInterpolationTable& interpTable
)
:
    basicInterpolationTable<Type>(interpTable)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationTables::linearInterpolationTable<Type>::
~linearInterpolationTable()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::List<Type> 
Foam::interpolationTables::linearInterpolationTable<Type>::interpolate
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

    // Extract upstream element
	const Tuple2<scalar, List<Type>>& us = this->values_[nextElement];
	if(us.first() == pTime)
	{
		return us.second();
	}

    // Extract downstream element
	const Tuple2<scalar, List<Type>>& ds = this->values_[nextElement-1];

	if (ds.second().size() != us.second().size())
	{
		FatalErrorInFunction
			<< "Elements " << nextElement-1 << " and "
            << nextElement << "of time series "
			<< " have different sizes ... Can't interpolate."
			<< exit(FatalError);
	}

    // Weight value
	scalar w  = (pTime-ds.first())/(us.first()-ds.first());

    // Return weighted list
	List<Type> result(us.second().size());
	forAll(result, i)
	{
		result[i] = (1-w)*ds.second()[i] + w*us.second()[i];
	}

	return result;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationTables::linearInterpolationTable<Type>&
Foam::interpolationTables::linearInterpolationTable<Type>::operator=
(
    const linearInterpolationTable& interpTable
)
{
    *this = Foam::basicInterpolationTable<Type>::operator=(interpTable);
    return *this;
}

// ************************************************************************* //
