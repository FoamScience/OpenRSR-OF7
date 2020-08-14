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

\*---------------------------------------------------------------------------*/

#include "linearInterpolationTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::linearInterpolationTable<Type>::linearInterpolationTable
(
    const List<Tuple2<scalar, List<Type>>>& values,
	const bool isPeriodic
)
:
    basicInterpolationTable<Type>(values, isPeriodic)
{
}

template<class Type>
Foam::linearInterpolationTable<Type>::linearInterpolationTable(const dictionary& dict)
:
    basicInterpolationTable<Type>(dict)
{
}


template<class Type>
Foam::linearInterpolationTable<Type>::linearInterpolationTable
(
     const linearInterpolationTable& interpTable
)
:
    basicInterpolationTable<Type>(interpTable)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::linearInterpolationTable<Type>::~linearInterpolationTable()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::List<Type> Foam::linearInterpolationTable<Type>::interpolate
(
    const scalar& time
) const
{
	if (this->values_.size() <= 1)
	{
		return this->values_[0].second();
	}
	
	scalar pTime = this->projectTime(time);
	int nextElement = this->lookup(pTime);

	const Tuple2<scalar, List<Type>>& us = this->values_[nextElement];

	if(us.first() == pTime)
	{
		return us.second();
	}

	const Tuple2<scalar, List<Type>>& ds = this->values_[nextElement-1];

	if (ds.second().size() != us.second().size())
	{
		FatalErrorInFunction
			<< "Elements " << nextElement-1 << " and "
            << nextElement << "of time series "
			<< " have different sizes ... Can't interpolate."
			<< exit(FatalError);
	}

	scalar w  = (pTime-ds.first())/(us.first()-ds.first());
	List<Type> result(us.second().size());
	forAll(result, i)
	{
		result[i] = (1-w)*ds.second()[i] + w*us.second()[i];
	}

	return result;
}


// ************************************************************************* //
