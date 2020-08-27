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
    Foam::steppedInterpolationTable

Description
    Time series with step-like interpolation.

SourceFiles
    steppedInterpolationTable.C

\*---------------------------------------------------------------------------*/

#ifndef steppedInterpolationTable_H
#define steppedInterpolationTable_H

#include "basicInterpolationTable.H"
#include "tableReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interpolationTables
{

/*---------------------------------------------------------------------------*\
                     Class steppedInterpolationTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class steppedInterpolationTable
:
    public basicInterpolationTable<Type>
{

public:

    //- Runtime type information
    TypeName("stepped");

    // Constructors

        //- Construct from components
        steppedInterpolationTable
        (
            const List<Tuple2<scalar, List<Type>>>& values,
			const bool isPeriodic
        );

        //- Construct from dictionary
        steppedInterpolationTable(const dictionary& dict);

        //- Construct copy
        steppedInterpolationTable(const steppedInterpolationTable& interpTable);

        //- Construct and return a copy
        virtual autoPtr<basicInterpolationTable<Type>> clone() const
        {
            return autoPtr<basicInterpolationTable<Type>>
            (
                new steppedInterpolationTable<Type>
                (
                    *this
                )
            );
        }

    //- Destructor
    virtual ~steppedInterpolationTable();

    // Member Functions

		//- The actual interpolation operation
		virtual List<Type> interpolate(const scalar& time) const;

    // Member Operators
        
        //- Disallow default bitwise assignment
        steppedInterpolationTable& operator=(const steppedInterpolationTable&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interpolationTables

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "steppedInterpolationTable.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
