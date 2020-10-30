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
    Foam::capPressModel

Description
    A template base class for rock-dependent capillary pressure models.

SourceFiles
    capPressModel.C

\*---------------------------------------------------------------------------*/

#ifndef capPressModel_H
#define capPressModel_H

#include "DimensionedField.H"
#include "phase.H"
#include "surfaceMesh.H"
#include "rock.H"
#include "fvMesh.H"
#include "HashTable.H"
#include "addToTemplatedRunTimeSelection.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Declare rock types
using iRock = rock<Isotropic>;
using dRock = rock<DiagAnisotropic>;
using aRock = rock<Anisotropic>;

/*---------------------------------------------------------------------------*\
                         Class capPressModel Declaration
\*---------------------------------------------------------------------------*/

template<class RockType, int nPhases>
class capPressModel
:
    public regIOobject
{
public:

    // Public Member Variables

        // The pattern controlling suggested model names
        static string nameRegExp;

    // Static Member Functions

        //- Parse Model name to (example: pcModel<phase1,phase2>)
        // 
        // - First: the parsed pure name (pcModel)
        // - Second: a list of phase names (phase1, phase2)
        static Tuple2<word,wordList> parseModelName(const word& modelName);

        //- Return Pc field name for given phase
        static word pcName(const word& phaseName)
        {
            return phaseName+".pc";
        }

        //- Return pc derivative field name for given phase wrt canonical phase
        static word dpcName
        (
            const word& phaseName, const word& canonicalPhase
        )
        {
            return phaseName+".dpcdS("+canonicalPhase+")";
        }

        //- Pick the first model of a list which calculates some pc for
        //  the given phase
        static const capPressModel& getPcModel
        (
            const word& phaseName,
            const fvMesh& mesh
        );

protected:

    // Protected Data

        //- Model name
        word name_;

        //- Model Dict
        dictionary pcDict_;

        //- Const access to the rock object
        const RockType& rock_;

        //- Target Phase Names for this Pc model
        wordList phaseNames_;

        //- Canonical phase names (Those with increasing saturation)
        wordList canonicalPhases_;

        //- Phase relative permeabilities and their derivatives
        //  in respect to the canonical phases.
        HashTable<volScalarField> pcTable_;

public:

    //- Runtime type information
    ClassName("capPressModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        capPressModel,
        dictionary,
        (
            const word& name,
            const dictionary& transportProperties,
            const RockType& rock
        ),
        (name, transportProperties, rock)
    );

    // Constructors

        //- Construct from components
        capPressModel
        (
            const word& name,
            const dictionary& transportProperties,
            const RockType& rock
        );

        //- Construct from copy
        capPressModel(const capPressModel&) = delete;

    //- Destructor
    virtual ~capPressModel();

    // Selectors

        //- Return a reference to the selected capPressModel
        static autoPtr<capPressModel> New
        (
            const word& name,
            const dictionary& transportProperties,
            const RockType& rock
        );

    // Member Functions

        //- Return Pc model name
        const word& name() const {
            return name_;
        }

        //- Return Pc Model dictionary
        const dictionary& dict() const {
            return pcDict_;
        }

        //- Return const-access to the rock
        const RockType& rock() const {
            return rock_;
        }

        //- Return phase names list
        const wordList& phases() const {
            return phaseNames_;
        }

        //- Return canonical phases
        const wordList& canonicalPhases() const {
            return canonicalPhases_;
        }

        //- Write Pc data to Ostream
        virtual bool writeData(Ostream&) const override
        {
            NotImplemented;
            return false;
        }

        //- Update Pc table
        virtual void correct() = 0;

    // Member Operators

        //- Write-Access to pc fields from table
        volScalarField& operator[]
        (
            const word& fieldName
        )
        {
            if (!pcTable_.found(fieldName))
            {
                FatalErrorInFunction
                    << "No field named " << fieldName
                    << "in capillary pressure table."
                    << endl << "Valid fields are "
                    << pcTable_.sortedToc();
            }
            return pcTable_[fieldName];
        }

        //- Const-access to pc fields from table
        const volScalarField& operator[]
        (
            const word& fieldName
        ) const
        {
            if (!pcTable_.found(fieldName))
            {
                FatalErrorInFunction
                    << "No field named " << fieldName
                    << "in capillary pressure table."
                    << endl << "Valid fields are "
                    << pcTable_.sortedToc();
            }
            return pcTable_[fieldName];
        }

        //- Disallow default bitwise assignment
        capPressModel& operator=(const capPressModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#define maketwoPhaseCapPressTemplatedModel\
(baseType,childType,TargName,ns,childTArg)                                     \
                                                                               \
    /* Define child type name */                                               \
    template<>                                                                 \
    const ::Foam::word ns::childType<childTArg >::                             \
        typeName(ns::childType<childTArg >::typeName_());                      \
                                                                               \
    /* Define child debug switch */                                            \
    template<>                                                                 \
    int ns::childType<childTArg >::debug                                       \
    (                                                                          \
        ::Foam::debug::debugSwitch(ns::childType<childTArg >::typeName_(), 0)  \
    );                                                                         \
    /* Register the debug switch */                                            \
    template<>                                                                 \
    const Foam::RegisterDebugSwitch<ns::childType<childTArg >>                 \
    Foam::RegisterDebugSwitch<ns::childType<childTArg >>                       \
        ::registerDebugSwitch(ns::childType<childTArg >::typeName_());         \
                                                                               \
    /* Add child to the base's run-time selection tabe */                      \
    namespace ns {                                                             \
    baseType<childTArg, 2>::                                                   \
    adddictionaryConstructorToTable<childType<childTArg >>                     \
    add##phaseType##TargName##argNames##ConstructorTo##baseType##Targ##Table_; \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "capPressModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
