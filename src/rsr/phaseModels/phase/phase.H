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

#include "error.H"
#include "surfaceFieldsFwd.H"
#include "volFields.H"
#include "regIOobject.H"
#include "fvcFlux.H"
#include "FVFModel.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcReconstruct.H"
#include "UniformityTypes.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class phase Declaration
\*---------------------------------------------------------------------------*/

// Supported mixture types
enum class mixtureType { singlePhase, multiPhase };

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

        //- Phase saturation in a mixture
        autoPtr<volScalarField> alphaPtr_;

        //- Standard phase density
        dimensionedScalar rhoSc_;

        //- Phase density
        CompressibilityType rho_;

        //- Phase flux
        autoPtr<surfaceScalarField> phiPtr_;

        //- Phase viscosity
        ViscosityType mu_;

        //- FVF Model
        autoPtr<FVFModel<CompressibilityType>> BModel_;

public:

    //- Runtime type information
    ClassName("phase");

    // Constructors

        // Construct from components
        phase
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const mixtureType& mT = mixtureType::multiPhase
        );

        // Construct from copy
        phase (const phase&);


    //- Destructor
    virtual ~phase();


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

        //- Return const-ref to density field
        const CompressibilityType& rho() const {
            return rho_;
        }

        //- Return const-ref to viscosity field
        const CompressibilityType& mu() const {
            return mu_;
        }

        //- Return const-access to FVF Model
        const FVFModel<CompressibilityType>& BModel()
        {
            return BModel();
        }

        //- Write phase data to Ostream
        virtual bool writeData(Ostream&) const override
        {
            NotImplemented;
            return false;
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

#ifdef NoRepository
    #include "phase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
