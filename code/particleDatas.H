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
    Foam::particleDatas

Description

SourceFiles
    particleDatasI.H
    particleDatas.C
    particleDatasIO.C

\*---------------------------------------------------------------------------*/

#ifndef PARTICLEDATAS_H
#define PARTICLEDATAS_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
// Forward declaration of friend functions and operators


/*---------------------------------------------------------------------------*\
                         Class CLASSNAME Declaration
\*---------------------------------------------------------------------------*/

class particleDatas

{
public:
    // Public data
	//QMOM data
// cell id
    label cid_;
// number of volume node
    const label Nv_;
// number of area node
//    const label Na_;
// number of fractal dimension node
    const label Nf_;
//number of moments
    //volume
    const label nMomentv_;
    //fractal dimension
    const label nMomentf_;
    //total
    const label nMoment_;
    //total node number
    const label nNode_;
// node id
    const label nodeID_;
// volume abscissa
    scalar volume_;
// area abscissa
    scalar area_;
// fractal dimension abcissa
    scalar Df_;
//asymptotic fractal dimension
    const scalar Dfasymp_;
// fractal prefactor
//    scalar k_;       ////どこに入れる？

// volume weight
    scalar weightv_;
// area weight
//    scalar weighta_;
// fractal dimension weight under each volume
    scalar weightf_;

    //Field data
// temperature
    scalar T_;
// density
    scalar rho_;
// viscosity
    scalar mu_;
// mean free path of the gas
    scalar lambdaG_;

    //particle data
// equivalent primary paricle diameter
    scalar Dp_;
// number of primary particles in agglomerate
    scalar Np_;
// particle collision diameter
    scalar Dc_;
// diffusion coefficient of each particles
    scalar D_;
// particle velocity
    scalar c_;
// Knudsen number
    scalar Kn_;
// mean free path of the particle
    scalar lambdaP_;
// radius of gyration
    scalar gy_;
// surface area of fully sintered particle
    scalar am_;
// grain bounary distance * grain boundary diffusion coefficient
    scalar deltaDb_;
// collision kernel against other particles
    scalarList betaC_;
// characteristic sintering time
    scalar tauS_;
// characteristic collision time against other particles
    scalarList tauC_;
// characteristic surface growth time
//    scalar tauSG_;
// nucleation rate
//    scalar nucRate_;
// nuclei size
//    scalar nucDiam_;
// SiO2 weight concentration
//    scalar SIO2G_;

// momentVF decrease against all particles by collision
//    scalarList collision_;

// momentVF decrease against all particles by collision
    scalarRectangularMatrix collisionVF_;
// momentVF decrease against all particles by collision
    scalarRectangularMatrix collisionA_;
    //        dimensionedScalar tauCMin_;


    // Private Member Functions

    //- Disallow default bitwise copy construct
    //        particleDatas(const particleDatas&);

    //- Disallow default bitwise assignment
    //        void operator=(const particleDatas&);


public:

    // Static data members

    static const particleConstants CONST_;
    static label nextNodeID;


    // Constructors

    //- Construct null
    particleDatas();

    //- Construct from nodeID
    //        particleDatas(const label& nodeID);

    //- Construct from components

    //- Construct from Istream
    //        particleDatas(Istream&);

    //- Construct as copy
    //        particleDatas(const particleDatas&);


    // Selectors

    //- Select null constructed
    //        static autoPtr<particleDatas> New();


    //- Destructor
    ~particleDatas();


    // Member Functions

    // Access
    label& cid();
    const scalar& collision(const label& iCol) const;
    // dimensionedScalar -> scalar
    scalar& volume();
    scalar& area();
    scalar& Df();
    scalar& weightv();
    scalar& weighta();
    scalar& weightf();
    const scalar& Dp() const;
    const scalar& Dc() const;
    const scalar& D() const;
    const scalar& tauS() const;
    const scalar& betaC(const label& iOpp) const;

    // Check

    // Edit

    // get weight and abcissa value
    void getWeightsAbscissas(
    		const PtrList<volScalarField>& weightv,
            const PtrList<volScalarField>& weightf,
            const PtrList<volScalarField>& abscissav,
            const PtrList<volScalarField>& abscissaf,
            const PtrList<volScalarField>& abscissaa
            //const PtrList<volScalarField>& Lac
            //const PtrList<volScalarField>& weighta,
    );

    // set each field value
    void setcid(
            const label& cid,
            const volScalarField& T,
            const volScalarField& mu,
            const volScalarField& lambdaG,
            const volScalarField& deltaDb,
           // const volScalarField& SIO2G
            const volScalarField& rho
//            const volScalarField& nucRate,
//            const volScalarField& nucDiam
    );

    //set particle data
    void set();

    //update collision kernel
    void updateBetaC(const List<particleDatas>& DatasList);

    //calculate tauC
    scalar tauCMin(const List<particleDatas>& DatasList);

    //calculate surface area variation by sintering
    void sinteringRestructuring(const scalar& deltaTPLoop);

    //calculate surface area & volume variation by surface growth
    //void surfaceGrowth(const scalar& deltaTPLoop);

    //calculate nucleation rate & nuclei size
    //void nucleation(const scalar& deltaTPLoop);
    // Write


    // Member Operators

    //        void operator=(const particleDatas&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

    //        friend Istream& operator>>(Istream&, particleDatas&);
    //        friend Ostream& operator<<(Ostream&, const particleDatas&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "particleDatasI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif //PARTICLEDATAS_H

// ************************************************************************* //
