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
    Foam::functionObjects::actualWellFlowrate

Description
    Calculates the actual (per-phase) flowrate into/out of user-specified wells
    directly from well matrices.

    This operation can be carried out at reservoir (default) or surface
    conditions.

    Example of function object specification:
    \verbatim
    actualWellFlowrate1
    {
        type        actualWellFlowrate;
        libs        ("librsrFunctionObjects.so");
        write       yes;
        log         yes;
        collated    yes;
        conditions  surface;
        wellModel   wModel;
        wells
        (
            defaultGrp
            {
                phases ( water oil );
                wellNames ( PROD0 INJ0 );
            }
            group1
            {
                phases ( water );
                wellNames (PROD2);
            }
        );
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: actualWellFlowrate  | yes         |
        write        | write min/max data to file |  no      | yes
        log          | write min/max data to standard output | no | no
        collated     | write each well/phase data to a seperate file | no | no
        conditions   | calculate flowrates at reservoir or surface conditions | no | "reservoir"
        wellModel    | well model name to target | no | "wModel"
        wells        | list of groupped well names with requested phases | yes | |
    \endtable

    Output data is written to the file \<timeDir\>/actualWellFlowrate.dat

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles

SourceFiles
    actualWellFlowrate.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_actualWellFlowrate_H
#define functionObjects_actualWellFlowrate_H

#include "fvMeshFunctionObject.H"
#include "wellModel.H"
#include "logFiles.H"
#include "vector.H"
#include "rock.H"
#include "entry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class actualWellFlowrate Declaration
\*---------------------------------------------------------------------------*/

class actualWellFlowrate
:
    public fvMeshFunctionObject,
    public logFiles
{
public:

    enum class modeType
    {
        reservoir,
        surface
    };

protected:

    // Protected data

        //- Mode type names for conditions
        static const NamedEnum<modeType, 2> modeTypeNames_;

        //- Mode for conditions at which calculations are carried out
        modeType mode_;

        //- Write each well per-phase rate to a sperate file (Default: false)
        bool collated_;

        //- Well Model Object name
        word wModelName_;

        //- Requested Rates table
        HashTable<scalar, word> rates_;

        //- Was well source for a phase was calculated already in this timeStep?
        HashTable<bool, word> sourceNotCalculated_;

        //- Calculated well sources
        HashTable<scalarList, word> sources_;

        //- Output file names
        wordList fileNames_;


    // Protected Member Functions

        //- Calculate actual flowrate
        template<class RockType, int nPhases>
        void calcActualRate
        (
            const word& rateName,
            const modeType& mode
        );

        //- Output file header information
        virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("actualWellFlowrate");


    // Constructors

        //- Construct from Time and dictionary
        actualWellFlowrate
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        actualWellFlowrate(const actualWellFlowrate&) = delete;


    //- Destructor
    virtual ~actualWellFlowrate();


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the actualWellFlowrate
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const actualWellFlowrate&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "actualWellFlowrateTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
