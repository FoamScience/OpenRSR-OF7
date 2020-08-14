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
    Foam::linearInterpolationTable

Description
    An abstract interpolation/look-up TimeSeries of scalar vs 
    \<ElementType\> values.
    The reference time scalar values must be monotonically increasing.

    If \c isPeriodic is enabled, the final time value is treated as 
	being equivalent to start time for the following periods.

    Read csv format:
    \verbatim
        readerType      csv;
        file            "$FOAM_CASE/constant/krVsSat.csv";
        hasHeaderLine   true;   // skip first line
        timeColumn      0;      // time is in column 0
        valueColumns    (1);    // value starts in column 1
    \endverbatim


Note
    - Empty TimeSeries results in a FatalError

SourceFiles
    linearInterpolationTable.C

\*---------------------------------------------------------------------------*/

#ifndef linearInterpolationTable_H
#define linearInterpolationTable_H

#include "basicInterpolationTable.H"
#include "tableReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class linearInterpolationTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class linearInterpolationTable
:
    public basicInterpolationTable<Type>
{
public:

    //- Runtime type information
    TypeName("linearInterpolationTable");


    // Constructors

        //- Disallow null constructor
        linearInterpolationTable() = delete;

        //- Construct from components
        linearInterpolationTable
        (
            const List<Tuple2<scalar, List<Type>>>& values,
			const bool isPeriodic
        );

        //- Construct by reading the fileName from dictionary
        //  and read the table from that file.
        explicit linearInterpolationTable(const dictionary& dict);

        //- Construct copy
        linearInterpolationTable(const linearInterpolationTable& interpTable);

	// Destructor
	virtual ~linearInterpolationTable();

    // Member Functions

		//- The actual interpolation operation
		virtual List<Type> interpolate(const scalar& time) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "linearInterpolationTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
