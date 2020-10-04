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

\*---------------------------------------------------------------------------*/

#ifndef wellI_H
#define wellI_H

#include "well.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RockType, int nPhases>
inline void Foam::well<RockType, nPhases>::cellsVolume()
{
    V_.value() = 0;
    forAll(wellSet_.toc(), ci)
    {
        V_ += rock_.mesh().V()[wellSet_.toc()[ci]];
    }
    //- V_ is now the well cells volume managed by this processor.
    //- Reduce and scatter V_ to send the true total value to all processors
    reduce(V_, sumOp<scalar>());
    Pstream::scatter(V_);
}

template<class RockType, int nPhases>
inline Foam::label Foam::well<RockType, nPhases>::operationSign() const
{
    return  operation_ == operationHandling::production ? 1 : -1;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
