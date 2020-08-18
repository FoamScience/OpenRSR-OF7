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

#include "sampleInterpolationTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationTables::sampleInterpolationTable<Type>::
sampleInterpolationTable
(
    const List<Tuple2<scalar, List<Type>>>& values,
	const bool isPeriodic
)
:
    basicInterpolationTable<Type>(values, isPeriodic)
{
}

template<class Type>
Foam::interpolationTables::sampleInterpolationTable<Type>::
sampleInterpolationTable
(
    const dictionary& dict
)
:
    basicInterpolationTable<Type>(dict)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationTables::sampleInterpolationTable<Type>::
~sampleInterpolationTable() {}


template<>
Foam::List<Foam::scalar>
Foam::interpolationTables::sampleInterpolationTable<Foam::scalar>::interpolate
(
    const scalar& time
) const
{                
    return
    {
        this->projectTime(time),
        static_cast<scalar>(this->lookup(time))
    };
}

// ************************************************************************* //
