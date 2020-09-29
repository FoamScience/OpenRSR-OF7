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
    A generic base class for fluid phases supporting incompressible and
    compressible (through FVFModel) setups in single and multi-phase 
    contexts.

    Compressibility handling is the responsibility of FVFModel class.
        
    Optimizing viscosity field size is subject to the following reasoning:
    - If a mu0 dimensionedScalar is present in phase dictionary, mu is
      initialized to have only one element (the value is the provided one).
    - If not, 0/phasename.mu must be read from case files and it assumes
      full mesh.

Note
    0/phasename.mu is always preferred over dimensionedScalar from phase dict.

SourceFiles
    phase.C

\*---------------------------------------------------------------------------*/

#ifndef phase_H
#define phase_H

#include "VolatileDimensionedField.H"
#include "phase.H"
#include "FVFModel.H"
#include "UniformityTypes.H"
#include "singleCellFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Supported mixture types
enum class mixtureType { singlePhase, multiPhase };

/*---------------------------------------------------------------------------*\
                            Class phase Declaration
\*---------------------------------------------------------------------------*/

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

        //- Phase saturation in a mixture
        autoPtr<volScalarField> alphaPtr_;

        //- Standard phase density
        dimensionedScalar rhoSc_;

        //- Initial density
        dimensionedScalar rho0_;

        //- Initial viscosity
        dimensionedScalar mu0_;

        //- Mesh for viscosity field
        autoPtr<singleCellFvMesh> muMesh_;

        //- Phase flux
        autoPtr<surfaceScalarField> phiPtr_;

        //- FVF Model
        autoPtr<FVFModel> BModel_;

        //- Phase density
        VolatileDimensionedField<scalar,volMesh> rho_;

        //- Phase viscosity
        VolatileDimensionedField<scalar,volMesh> mu_;

public:

    //- Runtime type information
    ClassName("phase");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        phase,
        dictionary,
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const mixtureType& mT
        ),
        (name, mesh, transportProperties, mT)
    );

    // Constructors

        // Construct from components
        phase
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const mixtureType& mT = mixtureType::multiPhase
        );

        // Construct from copy
        phase(const phase&);


    //- Destructor
    virtual ~phase();

    // Selectors

        //- Return a reference to the selected IncompressiblePhase
        static autoPtr<phase> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const mixtureType& mT = mixtureType::multiPhase
        );


    // Member Functions

        //- Return phase name
        const word& name() const {
            return name_;
        }

        //- Return phase dictionary
        const dictionary& dict() const {
            return phaseDict_;
        }

        //- Return const-access to phase mesh
        const fvMesh& mesh() const {
            return mesh_;
        }

        //- Return ref to velocity field (e.g. for boundary manipulation)
        volVectorField& U() {
            return U_;
        }

        //- Return const-ref to velocity field
        const volVectorField& U() const {
            return U_;
        }

        //- Return ref to saturation field (e.g. for boundary manipulation)
        volScalarField& alpha() {
            if (alphaPtr_.empty())
            {
                FatalErrorInFunction
                    << "Trying to access non-existent alpha field for phase "
                    << name_ << ". Is this a single phase simulation?" << endl
                    << exit(FatalError);
            }
            return alphaPtr_();
        }

        //- Return const-ref to velocity field
        const volScalarField& alpha() const {
            if (alphaPtr_.empty())
            {
                FatalErrorInFunction
                    << "Trying to access non-existent alpha field for phase "
                    << name_ << ". Is this a single phase simulation?" << endl
                    << exit(FatalError);
            }
            return alphaPtr_();
        }

        //- Return ref to phase flux
        surfaceScalarField& phi()
        {
            return phiPtr_();
        }

        //- Return cosnt-access to phase flux
        const surfaceScalarField& phi() const
        {
            return phiPtr_();
        }

        //- Return const-ref to standard density value
        const dimensionedScalar& rhoSc() const {
            return rhoSc_;
        }

        //- Write phase data to Ostream
        virtual bool writeData(Ostream&) const override
        {
            NotImplemented;
            return false;
        }

        //- Return const-ref to density field
        const VolatileDimensionedField<scalar, volMesh>& rho() const {
            return rho_;
        }

        //- Return const-ref to viscosity field
        const VolatileDimensionedField<scalar, volMesh>& mu() const {
            return mu_;
        }

        //- Return const-access to FVF Model
        const FVFModel& BModel() const
        {
            return BModel_();
        }

        //- Correct phase fields
        virtual void correct() = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        phase& operator=(const phase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "phase.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
