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
    Foam::phase

Description
    A template base class for fluid phases.

SourceFiles
    phase.C

\*---------------------------------------------------------------------------*/

#ifndef phase_H
#define phase_H

#include "volFields.H"
#include "regIOobject.H"
#include "fvcFlux.H"
#include "FVFModel.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcReconstruct.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class phase Declaration
\*---------------------------------------------------------------------------*/

template<class CompressibilityType, class ViscosityType>
class phase
:
    public regIOobject
{
protected:

    // Protected Data

        //- Phase name
        word name_;

        //- Phase Dict
        dictionary phaseDict_;

        //- Const access to the mesh
        const fvMesh& mesh_;

        //- Phase velocity
        volVectorField U_;

        //- Phase flux
        autoPtr<surfaceScalarField> phiPtr_;

        //- Phase saturation in a mixture
        volScalarField alpha_;

        //- Phase density
        CompressibilityType rho_;

        //- Phase viscosity
        ViscosityType mu_;

        //- FVF Model
        autoPtr<FVFModel> BModel_;

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
    TypeName("phase");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        phase,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    // Constructors

        //- Null constructor
        phase() = default;

        // Construct from components
        phase
        (
            const List<Tuple2<scalar, List<Type>>>& values,
            const bool isPeriodic
        );

        //- Construct from dictionary
        phase(const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<phase<Type>> clone() const = 0;

    // Selectors

        //- Return a reference to the selected phase
        static autoPtr<phase> New(const dictionary& spec);


    //- Destructor
    virtual ~phase();

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
        void operator=(const phase&);

        //- Return an element from values as const
        const Tuple2<scalar, List<Type>>& operator[](const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "phase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
