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

Class
    Foam::particleConstants

Description

SourceFiles
    particleConstantsI.H
    CLASSNAME.C
    CLASSNAMEIO.C

\*---------------------------------------------------------------------------*/

#ifndef PARTICLECONSTANTS_H
#define PARTICLECONSTANTS_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class particleConstants Declaration
\*---------------------------------------------------------------------------*/


struct particleConstants

{
    // Public data
	//nuclei's size is variable
    //const scalar DpNew;

    //const scalar V0;

    //const scalar A0;

	//volume nodes number
    const label Nv;

    //fractal dimension nodes number
    const label Nf;

    //fractal dimension of nuclei
    const scalar Df0;

    //fractal prefactor
    const scalar kf;

//    const label Na;

    //number of moments
    //volume
    const label nMomentv;
    //fractal dimension
    const label nMomentf;
    //total
    const label nMoment;
    //total node number
    const label nNode;

    //Boltzmann constant
    const scalar kb;

    //gas constant
    const scalar Rg;

	//density of SIO2
    const scalar rho_p;

	//molecular weight of SIO2
    const scalar M_p;

    //molecular volume of SiO2
    const scalar v_mol;

	//molecular surface area of SIO2
    const scalar a_mol;

	//molecular diameter of SIO2
    const scalar d_mol;

	//diameter of air molecular
    const scalar da;


    const scalar sigma;

    const scalar omega;

    const scalar Ea;

    const scalar deltaD0;

    const scalar pOne;

    const scalar Df_asymp;

    //surface tension of SIO2
    //const scalar st;

//    volScalarField deltaDb;

    // Static data members

    // Constructors

        //- Construct null
        particleConstants();


    // Selectors

        //- Select null constructed

    //- Destructor
    ~particleConstants();


    // Member Functions

        // Access

        // Check

        // Edit
//void setdeltaDb(const volScalarField& T);
        // Write


    // Member Operators


    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "particleConstantsI.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
