/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.orgfavorable texture is depicted by favorable texture is depicted by packaging of packaging of
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
    Foam::rock

Description
    A template base class for reservoir rocks, with pressure-dependent
    compressibility and anisotropic media support

SourceFiles
    rock.C

\*---------------------------------------------------------------------------*/

#ifndef rock_H
#define rock_H

#include "DimensionedField.H"
#include "GeometricField.H"
#include "regIOobject.H"
#include "UniformityTypes.H"
#include "IsotropyTypes.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class rock Declaration
\*---------------------------------------------------------------------------*/

template<class PermeabilityType, class CompressibilityType>
class rock
:
    public regIOobject
{
public:

    // Public Data Types

        // Permeability field component type
        using KcmptType = typename PermeabilityType::value_type;

        // Compressibility field component type
        using CcmptType = typename CompressibilityType::cmptType;

protected:

    // Protected Data

        //- Rock name
        word name_;

        //- Rock Dict
        dictionary rockDict_;

        //- Const access to the mesh
        const fvMesh& mesh_;

        //- Initial porosity
        dimensionedScalar poroInit_;

        //- Initial Permeability
        dimensioned<KcmptType> permInit_;

        //- Initial Compressibility
        dimensioned<CcmptType> compInit_;

        //- Rock porosity field
        volScalarField porosity_;

        //- Rock absolute-permeability field
        PermeabilityType K_;

        //- Rock Compressibility field
        CompressibilityType c_;

public:

    //- Runtime type information
    ClassName("rock");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        rock,
        dictionary,
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& rockProperties
        ),
        (name, mesh, rockProperties)
    );

    // Constructors

        // Construct from components
        rock
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& rockProperties
        );

        // Construct from copy
        rock (const rock&);


    //- Destructor
    virtual ~rock();

    // Selectors

        //- Return a reference to the selected DiagAnisoRock
        static autoPtr<rock> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& rockProperties
        );


    // Member Functions

        //- Return rock name
        const word& name() const {
            return name_;
        }

        //- Return rock dictionary
        const dictionary& dict() const {
            return rockDict_;
        }

        //- Return const-access to rock mesh
        const fvMesh& mesh() const {
            return mesh_;
        }

        //- Return const-ref to porosity field
        const volScalarField& porosity() const {
            return porosity_;
        }

        //- Return const-ref to absolute permeability field
        const PermeabilityType& K() const {
            return K_;
        }

        //- Return const-ref to compressibility field
        const CompressibilityType& Cf() {
            return c_;
        }

        //- Write rock data to Ostream
        virtual bool writeData(Ostream&) const override
        {
            NotImplemented;
            return false;
        }

        //- Correct rock fields
        virtual void correct() = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        rock& operator=(const rock&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rock.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
