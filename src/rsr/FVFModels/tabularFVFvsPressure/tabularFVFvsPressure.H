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
    Foam::FVFModels::tabularFVFvsPressure

Description
    Tabular description of the relationship 1|FVF = f(P). Typically used
    for undersaturated isotherm reservoirs.

SourceFiles
    tabularFVFvsPressure.C

\*---------------------------------------------------------------------------*/

#ifndef tabularFVFvsPressure_H
#define tabularFVFvsPressure_H

#include "FVFModel.H"
#include "basicInterpolationTable.H"
#include "fvPatchField.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FVFModels
{

/*---------------------------------------------------------------------------*\
                     Class tabularFVFvsPressure Declaration
\*---------------------------------------------------------------------------*/

class tabularFVFvsPressure
:
    public FVFModel<Compressible>
{
    // Private Data Members

        //- Pressure field name
        word pName_;

        //- Const access to pressure field
        const volScalarField& p_;

        //- Time series of (rFVF drFVFdP) vs Pressure
        autoPtr<basicInterpolationTable<scalar>> rFVFseries_;

public:

    //- Runtime type information
    TypeName("tabularFVFvsPressure");

    // Constructors

        //- Construct from components
        tabularFVFvsPressure
        (
            const word& name,
            const dictionary& phaseDict,
            const fvMesh& mesh
        );

    // Destructor
    virtual ~tabularFVFvsPressure() {}

    // Member Functions

        //- Correct the FVF
        virtual void correct() override;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FVFModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
