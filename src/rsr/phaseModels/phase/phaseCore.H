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
    Foam::phaseCore

Description
    Template-Invariant bits of phase objects

SourceFiles
    phaseCore.C

\*---------------------------------------------------------------------------*/

#ifndef phaseCore_H
#define phaseCore_H

#include "surfaceFieldsFwd.H"
#include "volFields.H"
#include "regIOobject.H"
#include "fvcFlux.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcReconstruct.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class phaseCore Declaration
\*---------------------------------------------------------------------------*/

//- Supported mixture types
enum class mixtureType { singlePhase, multiPhase };

class phaseCore
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

        //- Standard phaseCore density
        dimensionedScalar rhoSc_;

        //- Phase flux
        autoPtr<surfaceScalarField> phiPtr_;

public:

    // Constructors

        // Construct from components
        phaseCore
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const mixtureType& mT = mixtureType::multiPhase
        );

        // Construct from copy
        phaseCore(const phaseCore&);

    //- Destructor
    virtual ~phaseCore();

    // Member Functions

        //- Return phaseCore name
        const word& name() const {
            return name_;
        }

        //- Return phaseCore dictionary
        const dictionary& dict() const {
            return phaseDict_;
        }

        //- Return const-access to phaseCore mesh
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
                    << "Trying to access non-existent alpha field for phaseCore "
                    << name_ << ". Is this a single phaseCore simulation?" << endl
                    << exit(FatalError);
            }
            return alphaPtr_();
        }

        //- Return const-ref to velocity field
        const volScalarField& alpha() const {
            if (alphaPtr_.empty())
            {
                FatalErrorInFunction
                    << "Trying to access non-existent alpha field for phaseCore "
                    << name_ << ". Is this a single phaseCore simulation?" << endl
                    << exit(FatalError);
            }
            return alphaPtr_();
        }

        //- Return ref to phaseCore flux
        surfaceScalarField& phi()
        {
            return phiPtr_();
        }

        //- Return cosnt-access to phaseCore flux
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

    // Member Operators

        //- Disallow default bitwise assignment
        phaseCore& operator=(const phaseCore&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "phaseCore.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
