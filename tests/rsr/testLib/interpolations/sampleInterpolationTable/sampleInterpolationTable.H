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
    Foam::sampleInterpolationTable

Description


SourceFiles
    sampleInterpolationTable.C

\*---------------------------------------------------------------------------*/

#ifndef sampleInterpolationTable_H
#define sampleInterpolationTable_H

#include "basicInterpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interpolationTables
{

/*---------------------------------------------------------------------------*\
                   Class sampleInterpolationTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class sampleInterpolationTable
:
    public basicInterpolationTable<Type>
{
public:

    //- Runtime type information
    TypeName("sampleInterpolationTable");


    // Constructors

        //- Disallow null constructor
        sampleInterpolationTable() = delete;

        //- Construct from components
        explicit sampleInterpolationTable
        (
            const List<Tuple2<scalar, List<Type>>>& values,
			const bool isPeriodic = false
        );

        //- Construct by reading the fileName from dictionary
        //  and read the table from that file.
        explicit sampleInterpolationTable(const dictionary& dict);

        //- Construct and return a copy
        virtual autoPtr<basicInterpolationTable<Type>> clone() const
        {
            return autoPtr<basicInterpolationTable<Type>>
            (
                new sampleInterpolationTable<Type>
                (
                    *this
                )
            );
        }

	// Destructor
	virtual ~sampleInterpolationTable();

    // Member Functions

		//- The actual interpolation operation
		virtual List<Type> interpolate(const scalar& time) const
        {
            return {};
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interpolationTables

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampleInterpolationTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
