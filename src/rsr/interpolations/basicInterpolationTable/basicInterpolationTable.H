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
    Foam::basicInterpolationTable

Description
    An abstract interpolation/look-up time series of scalar vs 
    \<List<Type>\> values.
    The reference time scalar values must be monotonically increasing.

    If \c isPeriodic is enabled, the final time value is treated as 
	being equivalent to start time for the following periods.

    Interpolate "linearly":
    \verbatim
        interpolationType linear; // Default
        readerType      openFoam; // Default, only option for now  
        file            "$FOAM_CASE/constant/krVsSat.csv";
        periodic        false;
    \endverbatim

    Supported file readers:
        - OpenFOAM (default)

SourceFiles
    basicInterpolationTable.C

\*---------------------------------------------------------------------------*/

#ifndef basicInterpolationTable_H
#define basicInterpolationTable_H

#include "fileName.H"
#include "List.H"
#include "Tuple2.H"
#include "dictionary.H"
#include "tableReaders.H"
#include "fieldTypes.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class basicInterpolationTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class basicInterpolationTable
{
protected:

    // Protected Data

        //- (Unexpanded copy of) File name
        fileName fileName_;

        //- Values reader
        autoPtr<tableReader<List<Type>>> reader_;

        //- Values as a List
        List<Tuple2<scalar, List<Type>>> values_;

		//- Treat table as periodic?
		bool isPeriodic_;

		//- Time bounds
		scalar startTime_, endTime_;

    // Protected Member Functions

        //- Read the table of data from file into list
        void readTable();

		//- Project time value if isPeriodic_ is enabled
		scalar projectTime(scalar time) const;

		//- Get Index of the 1st value 
        //  whose time is greater or equal to requested time
		int lookup(const scalar& time) const;

public:

    //- Runtime type information
    TypeName("basicInterpolationTable");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicInterpolationTable,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    // Constructors

        // Construct from components
        basicInterpolationTable
        (
            const List<Tuple2<scalar, List<Type>>>& values,
            const bool isPeriodic
        );

        //- Construct from dictionary
        basicInterpolationTable(const dictionary& dict);

        //- Construct copy
        basicInterpolationTable(const basicInterpolationTable& interpTable);

        //- Construct and return a clone
        virtual autoPtr<basicInterpolationTable<Type>> clone() const = 0;

    // Selectors

        //- Return a reference to the selected basicInterpolationTable
        static autoPtr<basicInterpolationTable> New(const dictionary& spec);


    //- Destructor
    virtual ~basicInterpolationTable();

    // Member Functions

        //- Check that values are monotonically increasing
        //  else, FatalError-out
        void checkMonotonicity() const;

        //- Write to output stream
        void write(Ostream& os) const;

        //- Cosnt-access to underlying values list
        const List<Tuple2<scalar, List<Type>>> values() const {
            return values_;
        }

        //- Return time left until next value in list
        scalar deltaT(const scalar& time) const;

		//- The actual interpolation operation
		virtual List<Type> interpolate(const scalar& time) const = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        basicInterpolationTable& operator=(const basicInterpolationTable&);

        //- Return an element from values as const
        const Tuple2<scalar, List<Type>>& operator[](const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicInterpolationTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
