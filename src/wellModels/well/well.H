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
    Foam::well

Description
    A template base class for well objects which can belong to many well groups.

SourceFiles
    well.C

\*---------------------------------------------------------------------------*/

#ifndef well_H
#define well_H

#include "IOobject.H"
#include "error.H"
#include "phase.H"
#include "rock.H"
#include "cellSet.H"
#include "faceSet.H"
#include "wellSource.H"
#include "regIOobject.H"
#include "driveHandler.H"
#include "HashPtrTable.H"
#include "topoSetSource.H"
#include "sourceProperties.H"
#include "addToTemplatedRunTimeSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Declare rock types
using iRock = rock<Isotropic>;
using dRock = rock<DiagAnisotropic>;
using aRock = rock<Anisotropic>;

// A concrete regIOobject for well groups
struct wellGrpIO : public regIOobject
{
    //- Construct from components
    explicit wellGrpIO(const IOobject& io) : regIOobject(io) {}

    //- Should not be called to write any data to disk
    virtual bool writeData(Ostream&) const
    {
        NotImplemented;
        return false;
    }
};

/*---------------------------------------------------------------------------*\
                             Class well Declaration
\*---------------------------------------------------------------------------*/

template<class RockType, int nPhases>
class well
:
    public List<autoPtr<regIOobject>>
{

protected:

    // Protected Data

        //- Well name
        word name_;

        //- Well Dict
        dictionary wellDict_;

        //- Const access to the rock object
        const RockType& rock_;

        //- Well groups
        wordList wellGroups_;

        //- Well perforation intervals
        PtrList<topoSetSource> perfos_;

        //- Well drives
        PtrList<driveHandler<RockType, nPhases>> drives_;

        //- Well cells topoSet
        cellSet wellSet_;

        //- Well faces topoSet to hold faces internal to the well.
        faceSet faces_;

        //- Well properties for source calculations
        sourceProperties srcProps_;

        
    // Protected Member Functions

        //- Operation mode sign to use in calculations
        inline label operationSign() const;

public:

    //- Runtime type information
    ClassName("well");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        well,
        dictionary,
        (
            const word& name,
            const dictionary& wellDict,
            const RockType& rock,
            HashTable<autoPtr<wellSource<RockType, nPhases>>>& sources,
            HashPtrTable<fvScalarMatrix>& matTable
        ),
        (name, wellDict, rock, sources, matTable)
    );

    // Constructors

        //- Construct from components
        well
        (
            const word& name,
            const dictionary& wellDict,
            const RockType& rock,
            HashTable<autoPtr<wellSource<RockType, nPhases>>>& sources,
            HashPtrTable<fvScalarMatrix>& matTable
        );

        //- Construct from copy
        well(const well&) = delete;

    //- Destructor
    virtual ~well();

    // Selectors

        //- Return a reference to the selected well
        static autoPtr<well> New
        (
            const word& name,
            const dictionary& wellDict,
            const RockType& rock,
            HashTable<autoPtr<wellSource<RockType, nPhases>>>& sources,
            HashPtrTable<fvScalarMatrix>& matTable
        );

    // Member Functions

        //- Return well name
        const word& name() const {
            return name_;
        }

        //- Return well dictionary
        const dictionary& dict() const {
            return wellDict_;
        }

        //- Return const-access to the rock
        const RockType& rock() const {
            return rock_;
        }

        //- Return injection phase name
        const word& injPhase() const 
        {
            return srcProps_.injPhase();
        }

        //- Return operation mode
        const sourceProperties::operationHandling& operation() const 
        {
            return srcProps_.operation();
        }

        //- Return source-related properties
        const sourceProperties& srcProps() const
        {
            return srcProps_;
        }


        //- Return radius value
        const dimensionedScalar& radius() const
        {
            return srcProps_.radius();
        }

        //- Return skin value
        const scalar& skin() const
        {
            return srcProps_.skin();
        }

        //- Return well index
        const scalarList& wellIndex()
        {
            return srcProps_.wellIndex();
        }

        //- Return write-access to wellSet
        cellSet& wellSet() {
            return wellSet_;
        }

        //- Return per-processor cell IDs
        labelList cellIDs() const {
            return wellSet_.toc();
        }

        //- Return per-processor internal face IDs
        labelList faceIDs() const {
            return faces_.toc();
        }

        //- Return total volume of well cells
        const dimensionedScalar& totalVolume() const
        {
            return srcProps_.V();
        }

        //- Register well to groups from dictionary
        void registerToGroups();

        //- Read perforations
        void readPerforations();

        //- Read well drives
        void readImposedDrives
        (
            HashTable<autoPtr<wellSource<RockType, nPhases>>>& sources,
            HashPtrTable<fvScalarMatrix>& matTable
        );

        //- Update well sources
        virtual void correct() = 0;
        //void correct()
        //{
        //    forAll(drives_, di)
        //    {
        //        drives_[di].correct();
        //    }
        //};

    // Member Operators

        //- Disallow default bitwise assignment
        well& operator=(const well&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#define makeTemplatedWell\
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
    Foam::RegisterDebugSwitch<ns::childType<childTArg, nPhases>>               \
        ::registerDebugSwitch(ns::childType<childTArg, nPhases>::typeName_()); \
                                                                               \
    /* Add child to the base's run-time selection tabe */                      \
    namespace ns {                                                             \
    baseType<childTArg, nPhases>::                                             \
    adddictionaryConstructorToTable<childType<childTArg, nPhases>>             \
    add##phaseType##TargName##argNames##ConstructorTo##baseType##Targ##Table_; \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wellI.H"

#ifdef NoRepository
    #include "well.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
