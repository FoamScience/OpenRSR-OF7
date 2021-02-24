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
    Foam::impesControl

Description
    IMPES control class. Provides time-loop control methods which exit the
    simulation once convergence criteria have been reached with support for
    specialized IMPES settings. Example usage:

    \verbatim
    impesControl impes(mesh);

    while (impes.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // solve ...
    }
    \endverbatim

    - Approximates dSdt from alpha equation
    - Then calls a CFL calculation method

SourceFiles
    impesControl.C

\*---------------------------------------------------------------------------*/

#ifndef impesControl_H
#define impesControl_H

#include "fluidSolutionControl.H"
#include "singleRegionConvergenceControl.H"
#include "surfaceFieldsFwd.H"
#include "CFLMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class impesControl Declaration
\*---------------------------------------------------------------------------*/

template<class RockType, int nPhases>
class impesControl
:
    public fluidSolutionControl,
    public singleRegionConvergenceControl
{
protected:

    // Protected Data

        //- Maximal allowed change in canonical saturation between time levels
        scalar dSMax_;

        //- CFL Method to calculate CFL number
        autoPtr<CFLMethod<RockType,nPhases>> CFLMethod_;

public:

    //- Runtime type information
    ClassName("impesControl");

    // Constructors

        //- Construct from components
        impesControl
        (
            const word& name,
            const surfaceScalarField& phi,
            const wordList& phaseNames,
            const RockType& rock,
            const word& algorithmName="IMPES"
        );

    //- Destructor
    virtual ~impesControl();

    // Member Functions

        // IO

            //- Read controls
            virtual bool read();


        // Evolution

            //- Time run loop
            bool run(Time& time);

            //- Time loop loop
            bool loop(Time& time);

        // Getters

            //- Retun Max allowed saturation change
            scalar dSMax() const
            {
                return dSMax_;
            }

        // Time step calculation

            //- Time-step saturation equation from phase equation.
            //  Assume :
            //  ddt(coefft, phase.alpha) + div(phase.flux) == phase.src
            template<class CoeffFieldType>
            scalar deltaTFromAlphaEquation
            (
                const surfaceScalarField& phi,
                const scalarField& explicitSource,
                const CoeffFieldType& coefft 
            ) const;

            //- Return DeltaT from Calculated CFL number
            scalar deltaT
            (
                scalar maxCo,
                scalar dtForAlphaEqn,
                scalar maxDeltaT
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "impesControl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
