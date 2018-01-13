/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleConstants::particleConstants()
:
//DpNew(0.523e-9),
//V0(M_PI*DpNew*DpNew*DpNew/6.),
//A0(M_PI*DpNew*DpNew),
Nv(3),
Nf(3),
Df0(3),
kf(1),
//Na(1),
nMomentv(Nv*2),
nMomentf(Nf*2),
nMoment(Nv*(3*Nf+1)),//number of moments used in calculation
nNode(Nv*Nf),

kb(1.38e-23),
Rg(8.31),
rho_p(2.65e3),
M_p(60.08e-3),
d_mol(4.4e-10),
v_mol(M_PI*d_mol*d_mol*d_mol/6.),
a_mol(M_PI*d_mol*d_mol),
da(3.70e-10),
sigma(1.),
omega(4.49e-5),
Ea(4.1e5),
deltaD0(2.36e-12),
pOne(1.0),
Df_asymp(1.78)


{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleConstants::~particleConstants()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
/*void Foam::particleConstants::setdeltaDb(const volScalarField& T)
{
    deltaDb = deltaD0*exp(-Ea/(Rg*T));
}
*/
// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
