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
    Foam::driveHandler

Description
    A template base class for imposed drive handlers.

    Resposibilities:
    - Figure out the coefficients a, b and c in the relationship:
      $\f q = ap + b BHP + c\f$ for each phase (p is the canonical pressure)
      for a well *when* a parameter is imposed.
    - Take into account the well's operation mode

SourceFiles
    driveHandler.C

\*---------------------------------------------------------------------------*/

#ifndef driveHandler_H
#define driveHandler_H

#include "fvMatricesFwd.H"
#include "scalarList.H"
#include "wellSource.H"
#include "cellToFace.H"
#include "HashPtrTable.H"
#include "sourceProperties.H"
#include "basicInterpolationTable.H"
#include "addToTemplatedRunTimeSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class driveHandler Declaration
\*---------------------------------------------------------------------------*/

template<class RockType, int nPhases>
class driveHandler
{

protected:

    // Protected Data

        //- Drive name
        word name_;

        //- Drive dictionary
        dictionary driveDict_;

        //- Well source describer
        HashTable<autoPtr<wellSource<RockType, nPhases>>>& wellSources_;

        //- Access to Source-related part of the well
        sourceProperties& srcProps_;

        //- Access to per-phase source matrices
        HashPtrTable<fvScalarMatrix>& matrices_;

        //- Well cells IDs (global)
        const labelList& cells_;

        //- Time series of drive values vs simulation time
        autoPtr<basicInterpolationTable<scalar>> driveSeries_;

        //- Well coeffs q = coeff0 p + coeff1 BHP + coeff2
        List<scalarList> coeffs_;

public:

    //- Runtime type information
    ClassName("driveHandler");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        driveHandler,
        dictionary,
        (
            const word& name,
            const dictionary& driveDict,
            HashTable<autoPtr<wellSource<RockType, nPhases>>>& sources,
            sourceProperties& srcProps,
            HashPtrTable<fvScalarMatrix>& matrices
        ),
        (name, driveDict, sources, srcProps, matrices)
    );

    // Constructors

        //- Construct from components
        driveHandler
        (
            const word& name,
            const dictionary& driveDict,
            HashTable<autoPtr<wellSource<RockType, nPhases>>>& sources,
            sourceProperties& srcProps,
            HashPtrTable<fvScalarMatrix>& matrices
        );

        //- Construct from copy
        driveHandler(const driveHandler&) = delete;

    //- Destructor
    virtual ~driveHandler();

    // Selectors

        //- Return a reference to the selected driveHandler
        static autoPtr<driveHandler> New
        (
            const word& name,
            const dictionary& driveDict,
            HashTable<autoPtr<wellSource<RockType, nPhases>>>& sources,
            sourceProperties& srcProps,
            HashPtrTable<fvScalarMatrix>& matrices
        );

    // Member Functions

        //- Return drive model name
        const word& name() const {
            return name_;
        }

        //- Return drive dictionary
        const dictionary& dict() const {
            return driveDict_;
        }

        //- Return const-access to the well source describer
        const wellSource<RockType, nPhases>& sourceDescriber
        (
            const word& phaseName
        ) const {
            return wellSources_[phaseName];
        }

        //- Update drive sources
        virtual void correct() = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        driveHandler& operator=(const driveHandler&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#define makeWellDriveTemplatedModel\
(baseType,childType,TargName,ns,nPhases,childTArg)                             \
                                                                               \
    /* Define child type name */                                               \
    template<>                                                                 \
    const ::Foam::word ns::childType<childTArg, nPhases>::                     \
        typeName(ns::childType<childTArg, nPhases>::typeName_());              \
                                                                               \
    /* Define child debug switch */                                            \
    template<>                                                                 \
    int ns::childType<childTArg, nPhases>::debug                               \
    (                                                                          \
        ::Foam::debug::debugSwitch                                             \
        (ns::childType<childTArg, nPhases>::typeName_(), 0)                    \
    );                                                                         \
    /* Register the debug switch */                                            \
    template<>                                                                 \
    const Foam::RegisterDebugSwitch<ns::childType<childTArg, nPhases>>         \
    Foam::RegisterDebugSwitch<ns::childType<childTArg,nPhases>>                \
        ::registerDebugSwitch(ns::childType<childTArg,nPhases>::typeName_());  \
                                                                               \
    /* Add child to the base's run-time selection tabe */                      \
    namespace ns {                                                             \
    baseType<childTArg, nPhases>::                                             \
    adddictionaryConstructorToTable<childType<childTArg, nPhases>>             \
    add##phaseType##TargName##argNames##ConstructorTo##baseType##Targ##Table_; \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "driveHandler.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
