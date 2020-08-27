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

#include "basicInterpolationTable.H"
#include "fileName.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::basicInterpolationTable<Type>>
Foam::basicInterpolationTable<Type>::New
(
    const dictionary& spec
)
{
    const word interpolationType = spec.lookupOrDefault<word>
    (
        "interpolationType",
        "linear"
    );

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(interpolationType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown interpolation type " << interpolationType
            << nl << nl
            << "Valid interpolation types : " << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<basicInterpolationTable<Type>>(cstrIter()(spec));
}

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::basicInterpolationTable<Type>::readTable()
{
    // Preserve the original (unexpanded) fileName to avoid absolute paths
    // appearing in Error reports
    fileName fName = fileName_;
    fName.expand();
    // Read data from file
    reader_()(fName, values_);

    if (values_.empty())
    {
        FatalErrorInFunction
            << "Table read from " << fName << " is empty" << nl
            << exit(FatalError);
    }
    // Check values monotonicity
    checkMonotonicity();
}


template<class Type>
Foam::scalar Foam::basicInterpolationTable<Type>::projectTime
(
	scalar time
) const
{
	if (time > endTime_ and !isPeriodic_)
	{
		FatalErrorInFunction
			<< "Out of time range: Got " << time
			<< " but time in time series ends at " << endTime_
			<< exit(FatalError);
	}
	while(time > endTime_)
	{
		time -= endTime_;
	}
	return time;
}


template<class Type>
int Foam::basicInterpolationTable<Type>::lookup
(
	const scalar& time
) const
{
	if (time < startTime_)
	{
		FatalErrorInFunction
			<< "Out of time range: Got " << time
			<< " but time in timeseries starts at " << startTime_
			<< exit(FatalError);
	}

	scalar pTime = projectTime(time);

	// Find the first element that has time is greater than pTime
	auto it = std::upper_bound(
		values_.begin(),
		values_.end(),
		pTime,
		[](scalar t, const Tuple2<scalar, List<Type>>& e1) { 
			return t <= e1.first();
		}
	);

	return std::distance(values_.begin(), it);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::basicInterpolationTable<Type>::basicInterpolationTable
(
    const List<Tuple2<scalar, List<Type>>>& values,
	const bool isPeriodic
)
:
    fileName_("fileNameIsUndefined"),
    reader_(),
    values_(values),
	isPeriodic_(isPeriodic),
	startTime_(values[0].first()),
	endTime_(values[values.size()-1].first())
{
	checkMonotonicity();
}


template<class Type>
Foam::basicInterpolationTable<Type>::basicInterpolationTable
(
    const dictionary& dict
)
:
    fileName_(dict.lookup("file")),
    reader_(tableReader<List<Type>>::New(dict)),
    values_(),
	isPeriodic_(dict.lookupOrDefault<bool>("periodic", false)),
	startTime_(0),
	endTime_(0)
{
    readTable();
    startTime_ = values_[0].first();
    endTime_ = values_[values_.size()-1].first();
}


template<class Type>
Foam::basicInterpolationTable<Type>::basicInterpolationTable
(
    const basicInterpolationTable& interpTable
)
:
    fileName_(interpTable.fileName_),
    reader_(interpTable.reader_->clone()),
    values_(interpTable.values_),
    isPeriodic_(interpTable.isPeriodic_),
    startTime_(interpTable.startTime_),
    endTime_(interpTable.endTime_)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::basicInterpolationTable<Type>::~basicInterpolationTable()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::basicInterpolationTable<Type>::checkMonotonicity() const
{
	auto it = std::adjacent_find(
		values_.begin(),
		values_.end(),
		[] (
            const Tuple2<scalar, List<Type>>& e1, 
            const Tuple2<scalar, List<Type>>& e2
        ) { 
			return e1.first() > e2.first();
		}
	);
    if (it != values_.end()) {
        FatalErrorInFunction
			<< "Time scale is not monotonically increasing."
			<< exit(FatalError);
    }
}


template<class Type>
void Foam::basicInterpolationTable<Type>::write(Ostream& os) const
{
	NotImplemented;
}

template<class Type>
Foam::scalar Foam::basicInterpolationTable<Type>::deltaT
(
    const scalar& time
) const
{
	scalar pTime = projectTime(time);

	// Get upper index of time in list
	int uInd = lookup(pTime);

	// If at the end, return 0
	if (uInd >= values_.size()-1)
	{
		return 0;
	}

	return values_[uInd+1].first() - pTime;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::basicInterpolationTable<Type>&
Foam::basicInterpolationTable<Type>::operator=
(
    const basicInterpolationTable& interpTable
)
{
    fileName_ = interpTable.fileName_;
    reader_ = interpTable.reader_;
    values_ = interpTable.values_;
	isPeriodic_ = interpTable.isPeriodic_;
	startTime_ = interpTable.startTime_;
	endTime_ = interpTable.endTime_;
}

template<class Type>
const Foam::Tuple2<Foam::scalar, Foam::List<Type>>&
Foam::basicInterpolationTable<Type>::operator[](const label i) const
{
    return values_[i];
}

// ************************************************************************* //
