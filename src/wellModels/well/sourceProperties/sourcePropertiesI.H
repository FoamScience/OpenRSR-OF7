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

#ifndef sourcePropertiesI_H
#define sourcePropertiesI_H

#include "sourceProperties.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline void Foam::sourceProperties::cellsVolume()
{
    V_.value() = 0;
    forAll(cells_, ci)
    {
        V_.value() += mesh_.V()[cells_[ci]];
    }
    //- V_ is now the well cells volume managed by this processor.
    //- Reduce and scatter V_ to send the true total value to all processors
    reduce(V_.value(), sumOp<scalar>());
    Pstream::scatter(V_.value());
}

inline Tuple2<label, label> Foam::sourceProperties::gLowerCell() const
{
    // Find lowest cell in this process
    label fCell = cells_[0];
    forAll(cells_, ci)
    {
        scalar cg = (mesh_.C()[ci] && g_).value();
        scalar base = (mesh_.C()[fCell] && g_).value();
        fCell = (cg >= base) ? cells_[ci] : fCell;
    }
    // Find lowest cell for all processes
    Tuple2<label,label> res;
    if (Pstream::parRun())
    {
        List<Tuple2<label,vector>> candidates(Pstream::nProcs());
        candidates[Pstream::myProcNo()] = 
            Tuple2<label,vector>(fCell, mesh_.C()[fCell]);
        if (Pstream::master())
        {
            Pstream::gatherList(candidates);
            vector gfPos = candidates[0].second();
            res.second() = 0;
            forAll(candidates, ci)
            {
                scalar cg = (candidates[ci].second() && g_).value();
                scalar base = (gfPos && g_).value();
                if (cg >= base)
                {
                    res.first() = candidates[ci].first();
                    gfPos = candidates[ci].second();
                    res.second() = ci;
                } 
            }
            Pstream::scatter(res);
        }
    } else {
        res.first() = fCell;
        res.second()= 0;
    }
    return res;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
