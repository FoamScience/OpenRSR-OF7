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

Namespace
    Foam::FVFModels

Description
    A namespace for Formation Volume Factor models.

Class
    Foam::FVFModel

Description
    A generic abstract class for Formation Volume Factor Calculation.

SourceFiles
    FVFModel.C

\*---------------------------------------------------------------------------*/

#ifndef FVFModel_H
#define FVFModel_H

#include "DimensionedField.H"
#include "volFieldsFwd.H"
#include "volFields.H"
#include "singleCellFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FVFModel Declaration
\*---------------------------------------------------------------------------*/

class FVFModel
{

protected:

    // Protected data

        //- Model Name
        word name_;

        //- Hosting dict
        dictionary phaseDict_;

        //- Const ref to mesh
        const fvMesh& mesh_;

        //- 1|FVF in the incompressible case
        dimensionedScalar incrFVF_;

        //- A One-cell mesh for optimized incompressible phase operations
        autoPtr<singleCellFvMesh> oneCellMesh_;

        //- 1|FVF field
        DimensionedField<scalar, volMesh> rFVF_;

        //- d(1|FVF)/dP field
        DimensionedField<scalar, volMesh> drFVFdP_;


public:

    //- Runtime type information
    TypeName("FVFModel");


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        FVFModel,
        dictionary,
        (
            const word& name,
            const dictionary& phaseDict,
            const fvMesh& mesh
        ),
        (name, phaseDict, mesh)
    );

    // Selectors

        //- Return a reference to the selected FVFModel model
        static autoPtr<FVFModel> New
        (
            const word& name,
            const dictionary& phaseDict,
            const fvMesh& mesh
        );


    // Constructors

        //- Construct from components
        FVFModel
        (
            const word& name,
            const dictionary& phaseDict,
            const fvMesh& mesh
        );

        //- Construct from copy
        FVFModel(const FVFModel&);

    //- Destructor
    virtual ~FVFModel() {}

    // Member Functions

        //- Return model name
        const word& name() const
        {
            return name_;
        }

        //- Return model hosting dictionary
        const dictionary& phaseDict() const
        {
            return phaseDict_;
        }

        //- Return registration mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return employed mesh
        const fvMesh& fieldsMesh() const
        {
            if (rFVF_.size() == 1)
            {
                return oneCellMesh_();
            }
            return mesh_;
        }

        //- Return const ref to 1|FVF
        const DimensionedField<scalar,volMesh>& rFVF() const
        {
            return rFVF_;
        }

        //- Return const ref d(1|FVF)/dP
        const DimensionedField<scalar,volMesh>& drFVFdP() const
        {
            return drFVFdP_;
        }

        //- Correct the FVF fields
        virtual void correct() = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        FVFModel& operator=(const FVFModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

//#ifdef NoRepository
//    #include "FVFModel.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
