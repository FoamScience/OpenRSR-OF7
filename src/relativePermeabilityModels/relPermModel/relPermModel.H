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
    Foam::relPermModel

Description
    A template base class for rock-dependent relative permeability models.

SourceFiles
    relPermModel.C

\*---------------------------------------------------------------------------*/

#ifndef relPermModel_H
#define relPermModel_H

#include "DimensionedField.H"
#include "phase.H"
#include "surfaceMesh.H"
#include "rock.H"
#include "fvMesh.H"
#include "HashTable.H"
#include "singleCellFvMesh.H"
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
                         Class relPermModel Declaration
\*---------------------------------------------------------------------------*/

template<class RockType, int nPhases>
class relPermModel
:
    public regIOobject
{
public:

    // Public Member Variables

        // The pattern controlling suggested model names
        static string nameRegExp;

    // Static Member Functions

        //- Parse Model name to (example: krModel<phase1, phase2>)
        // 
        // - First: the parsed pure name (krModel)
        // - Second: a list of phase names (phase1, phase2)
        static Tuple2<word,wordList> parseModelName(const word& modelName);

        //- Return Kr field name for given phase
        static word krName(const word& phaseName)
        {
            return phaseName+".kr";
        }

        //- Return kr derivative field name for given phase wrt canonical phase
        static word dkrName
        (
            const word& phaseName, const word& canonicalPhase
        )
        {
            return phaseName+".dkrdS("+canonicalPhase+")";
        }

protected:

    // Protected Data

        //- Model name
        word name_;

        //- Model Dict
        dictionary krDict_;

        //- Const access to the rock object
        const RockType& rock_;

        //- Target Phase Names for this Kr model
        wordList phaseNames_;

        //- Canonical phase names (Those with increasing saturation)
        wordList canonicalPhases_;

        //- Phase relative permeabilities and their derivatives
        //  in respect to the canonical phases.
        HashTable<DimensionedField<scalar, volMesh>> krTable_;

    // Protected Member Members
        
public:

    //- Runtime type information
    ClassName("relPermModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        relPermModel,
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
        relPermModel
        (
            const word& name,
            const dictionary& transportProperties,
            const RockType& rock
        );

        //- Construct from copy
        relPermModel(const relPermModel&) = delete;

    //- Destructor
    virtual ~relPermModel();

    // Selectors

        //- Return a reference to the selected relPermModel
        static autoPtr<relPermModel> New
        (
            const word& name,
            const dictionary& transportProperties,
            const RockType& rock
        );

    // Member Functions

        //- Return Kr model name
        const word& name() const {
            return name_;
        }

        //- Return Kr Model dictionary
        const dictionary& dict() const {
            return krDict_;
        }

        //- Return const-access to the rock
        const RockType& rock() const {
            return rock_;
        }

        //- Return canonical phases
        const wordList& canonicalPhases() const {
            return canonicalPhases_;
        }

        //- Write Kr data to Ostream
        virtual bool writeData(Ostream&) const override
        {
            NotImplemented;
            return false;
        }

        //- Update Kr table
        virtual void correct() = 0;

    // Member Operators

        //- Write-Access to kr fields from table
        DimensionedField<scalar, volMesh>& operator[]
        (
            const word& fieldName
        )
        {
            if (!krTable_.found(fieldName))
            {
                FatalErrorInFunction
                    << "No field named " << fieldName
                    << "in relative permeability table."
                    << endl << "Valid fields are "
                    << krTable_.sortedToc();
            }
            return krTable_[fieldName];
        }

        //- Const-access to kr fields from table
        const DimensionedField<scalar, volMesh>& operator[]
        (
            const word& fieldName
        ) const
        {
            if (!krTable_.found(fieldName))
            {
                FatalErrorInFunction
                    << "No field named " << fieldName
                    << "in relative permeability table."
                    << endl << "Valid fields are "
                    << krTable_.sortedToc();
            }
            return krTable_[fieldName];
        }

        //- Disallow default bitwise assignment
        relPermModel& operator=(const relPermModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#define maketwoPhaseRelPermTemplatedModel\
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
    #include "relPermModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
