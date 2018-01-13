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
const Foam::particleConstants Foam::particleDatas::CONST_;
Foam::label Foam::particleDatas::nextNodeID(0);
// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
/*
Foam::particleDatas::particleDatas()
:
    baseClassName(),
    data_()
{}
*/

Foam::particleDatas::particleDatas()
:
    cid_(0),
    Nv_(CONST_.Nv),
    Nf_(CONST_.Nf),
//    Na_(CONST_.Na),
    nodeID_(nextNodeID),
    volume_(CONST_.v_mol), //particle volume
    area_(CONST_.a_mol), //nucleiから持ってくる？
    Df_(CONST_.Df0),
    nMomentv_(CONST_.nMomentv),
    nMomentf_(CONST_.nMomentf),
    nMoment_(CONST_.nMoment),
    nNode_(CONST_.nNode),
    Dfasymp_(CONST_.Df_asymp),
//    k_(1),
    weightv_(0.),
    weightf_(0.),
//    weighta_(0.0),
    T_(298.),
    rho_(0.),
    mu_(0.),
    lambdaG_(0.),
    Dp_(CONST_.d_mol),
    Np_(1.),
    Dc_(CONST_.d_mol),
    D_(0.),
    c_(0.),
    Kn_(1.0),
    lambdaP_(0.0),
    gy_(0.0),
    am_(CONST_.a_mol),
    deltaDb_(1.0),
    betaC_(CONST_.nNode, 0.),
    tauS_(1.),
    tauC_(CONST_.nNode, 1.),
    //tauSG_(1.),

    //nucRate_(0),
	//nucDiam_(0),
	//SIO2G_(0),

    //collision_(CONST_.nNode,0.)
	collisionVF_(CONST_.nMomentv, CONST_.nMomentf, 0),
	collisionA_(CONST_.Nv, CONST_.Nf, 0)
    //検討中


{
	//increase(+1) when constructor called
++nextNodeID;

}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
/*
Foam::autoPtr<Foam::particleDatas>
Foam::particleDatas::New()
{
    return autoPtr<particleDatas>(new particleDatas);
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleDatas::~particleDatas()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
/*
label& Foam::particleDatas::cid()
{
    return cid_;
}

const scalar& Foam::particleDatas::collision(const label& iCol) const
{
    return collision_[iCol];
}

scalar& Foam::particleDatas::area()
{
    return area_.value();
}

scalar& Foam::particleDatas::volume()
{
    return volume_.value();
}

scalar& Foam::particleDatas::weightv()
{
    return weightv_.value();
}

scalar& Foam::particleDatas::weighta()
{
    return weighta_.value();
}

const scalar& Foam::particleDatas::Dp() const
{
    return Dp_.value();
}

const scalar& Foam::particleDatas::Dc() const
{
    return Dc_.value();
}

const scalar& Foam::particleDatas::D() const
{
    return D_.value();
}

const scalar& Foam::particleDatas::tau() const
{
    return tau_.value();
}

const scalar& Foam::particleDatas::betaC(const label& iOpp) const
{
    return betaC_[iOpp].value();
}
*/

void Foam::particleDatas::getWeightsAbscissas(  //OK
        const PtrList<volScalarField>& weightv,
        const PtrList<volScalarField>& weightf,
//        const PtrList<volScalarField>& weighta,
        const PtrList<volScalarField>& abscissav,
        const PtrList<volScalarField>& abscissaf,
        const PtrList<volScalarField>& abscissaa
//        const PtrList<volScalarField>& Lac
)
{
    weightv_ = weightv[nodeID_/CONST_.Nf][cid_];
    weightf_ = weightf[nodeID_][cid_];
//    weighta_ = weighta[nodeID_][cid_];
    volume_ = abscissav[nodeID_/CONST_.Nf][cid_];
    Df_ = abscissaf[nodeID_][cid_];
    area_ = abscissaa[nodeID_][cid_];
//    k_ = Lac[nodeID_][cid_];
}


void Foam::particleDatas::setcid( //OK
		const label& cid,
        const volScalarField& T,
        const volScalarField& mu,
        const volScalarField& lambdaG,
        const volScalarField& deltaDb,
       // const volScalarField& SIO2G,
        const volScalarField& rho
//        const volScalarField& nucRate,
//        const volScalarField& nucDiam
)
{
    cid_ = cid;
    T_ = T[cid];
    mu_ = mu[cid];
    deltaDb_ = deltaDb[cid];
    lambdaG_ = lambdaG[cid];
   // SIO2G_ = SIO2G[cid]; //SiO2 mass fraction
    rho_ = rho[cid];
//    nucRate_ = nucRate[cid];
//    nucDiam_ = nucDiam[cid];
}


void Foam::particleDatas::set()  //OK
{
	//surface area set
    scalar area = area_;
    if(area < 1.e-30)
    {
    	area = CONST_.a_mol;
        //Info << "smallest cid" << cid_ << endl;
    }

    //volume set
    scalar volume = volume_;
    if(volume < 1.e-35)
    {
        volume = CONST_.v_mol;
        //Info << "SMALLEST cid" << cid_ << endl;
    }

    //fractal dimension set
    scalar Df = Df_;
    if(Df < 1.)
    {
    	Df = 1.;
    	//Info << "over agg cid" << cid_ << endl;
    }
    else if(Df > CONST_.Df0)
    {
     	Df = CONST_.Df0;
        //Info << "OVER SPH cid" << cid_ << endl;
    }
    ////////
    area_ = area;
    volume_ = volume;
    Df_ = Df;
    //////////
    Dp_ = 6.*volume/area;
    Info << "OK" <<endl;
    Np_ = area*area*area/(36.*M_PI*volume*volume);
    Info << "OK" <<endl;
    Dc_ = Dp_*Foam::pow(Np_/CONST_.kf , 1/Df);
    Info << "OK" <<endl;
    c_ = sqrt(8.*CONST_.kb*T_/(M_PI*CONST_.rho_p*volume));
    Info << "OK" <<endl;
    Kn_ = 2.*lambdaG_/Dc_;
    Info << "OK" <<endl;
    D_ = CONST_.kb*T_*(5.+4.*Kn_+6*Kn_*Kn_+18.*Kn_*Kn_*Kn_)
                /(3.*M_PI*mu_*Dc_*(5.-Kn_+(8.+M_PI)*Kn_*Kn_));

}

void Foam::particleDatas::updateBetaC(const List<particleDatas>& DatasList)  //OK
{
	//initialization
	for(label icol = 0; icol < nMomentv_; ++icol)
	{
		for(label irow = 0; irow < nMomentf_; ++irow)
		{
			collisionVF_[irow][icol] = 0.;
		}
	}
	for(label icol = 0; icol < Nf_; ++icol)
	{
		for(label irow = 0; irow < Nv_; ++irow)
		{
			collisionA_[irow][icol] = 0.;
		}
	}


	//characteristic length
	lambdaP_ = 8*D_/(M_PI*c_);
	//transition parameter
	gy_ = (Foam::pow((Dc_+lambdaP_),3)-Foam::pow((Dc_*Dc_+lambdaP_*lambdaP_),1.5))/(3*Dc_*lambdaP_)-Dc_;

	//update betaC
	for(label i=0; i<nNode_; ++i)
	{
		if(i<nodeID_)
		{
			betaC_[i] = DatasList[i].betaC_[nodeID_];
		}
		else
		{
			//sum of diffucivity
			scalar DD = DatasList[i].D_+D_;
			//sum of collision diameter
			scalar DcDc = DatasList[i].Dc_+Dc_;
			//components of BetaC
			scalar gyAve2 = 2.*Foam::sqrt(DatasList[i].gy_*DatasList[i].gy_+gy_*gy_);
			scalar cAve = sqrt(DatasList[i].c_*DatasList[i].c_+c_*c_);
			//collision frecency(Seinfeld 1986)
			betaC_[i] = 2.*M_PI*DD*DcDc/(DcDc/(DcDc+gyAve2)+8.*DD/(cAve*DcDc));
		}

		//ここまでOK


		//compute characteristic collision time
		tauC_[i] = 1./(betaC_[i]*(DatasList[i].weightv_*DatasList[i].weightf_+CONST_.pOne));
	}

	//construct matrix to compute moment change
	scalarRectangularMatrix volMatrixBirth(nMomentv_, Nv_, 1);
	scalarRectangularMatrix volMatrixDeath(nMomentv_, Nv_, 1);
	scalarRectangularMatrix dfMatrixBirth(Nf_, Nv_, 0);
	scalarRectangularMatrix dfMatrixDeath(Nf_, Nv_, 0);
	scalarRectangularMatrix dfMatrixBirthForUpdate(Nf_, Nv_, 0);
	scalarList volOwnNode(nMomentv_, 1.);
	scalarList dfOwnNode(nMomentf_, 1.);


	//set data to volMatrixBirth & Death
	for(label icol=0; icol<Nv_; ++icol)
	{
		for(label irow=1; irow<nMomentv_; ++irow)
		{
			volMatrixDeath[irow][icol] = volMatrixDeath[irow-1][icol]*DatasList[Nf_*icol].volume_;
			volMatrixBirth[irow][icol] = volMatrixBirth[irow-1][icol]*(volume_+DatasList[Nf_*icol].volume_);
		}
	}
	//set data to dfMatrixBirth & Death
	for(label icol=0; icol<Nv_; ++icol)
	{
		for(label irow=0; irow<Nf_; ++irow)
		{
			label id = Nf_*icol+irow;
			dfMatrixDeath[irow][icol] = DatasList[id].Df_;

			//calculation new fractal dimension
			scalar DsumPow3 = (Np_*Dp_*Dp_*Dp_+DatasList[id].Np_*DatasList[id].Dp_*DatasList[id].Dp_*DatasList[id].Dp_)
									/(Np_+DatasList[id].Np_);
			scalar Dsum = pow(DsumPow3, 0.333);
			scalar Numerator = (Np_*Dp_*Dp_*Dp_+DatasList[id].Np_*DatasList[id].Dp_*DatasList[id].Dp_*DatasList[id].Dp_)/DsumPow3;
			scalar Denominator = (pow(Dp_, Dfasymp_)*pow(Np_, Dfasymp_/Df_) + pow(DatasList[id].Dp_, Dfasymp_)*pow(DatasList[id].Np_, Dfasymp_/DatasList[id].Df_))
										/pow(Dsum, Dfasymp_);
			dfMatrixBirth[irow][icol] = Dfasymp_*log(Numerator)/stabilise(log(Denominator),SMALL);
		}
	}
	dfMatrixBirthForUpdate = dfMatrixBirth;

	//set data to vol & Df of own node
	for(label k=1; k<nMomentv_; ++k)
	{
		volOwnNode[k] = volOwnNode[k-1]*volume_;
	}
	for(label l=1; l<nMomentf_; ++l)
	{
		dfOwnNode[l] = volOwnNode[l-1]*Df_;
	}

	//constant
	const scalar constantOfCollision = 0.5*weightv_*weightf_;

	//compute collision matrix
	for(label l=0; l<nMomentf_; ++l)
	{
		if(l==0) //volume only weighted momnet Mk00
		{
			for(label k=0; k<nMomentv_; ++k)
			{
				for(label id=0; id<nNode_; ++id)
				{
					scalar weight = betaC_[id]*DatasList[id].weightv_*DatasList[id].weightf_;
					scalar inc = weight*(volMatrixBirth[k][id/Nv_]);
					scalar dec1 = weight*volOwnNode[k];
					scalar dec2 = weight*volMatrixDeath[k][id/Nv_];
					collisionVF_[k][l] +=inc-dec1-dec2;

					if(k<Nv_)
					{
						collisionA_[k][l] += inc*(area_+DatasList[id].area_)-dec1*area_-dec2*DatasList[id].area_;
					}
				}
			}
		}
		else // volume & fractal dimension weighted
		{
			for(label k=0; k<Nv_; ++k)
			{
				for(label id=0; id<nNode_; ++id)
				{
					scalar weight = betaC_[id]*DatasList[id].weightv_*DatasList[id].weightf_;
					scalar inc = weight*(volMatrixBirth[k][id/Nv_]*dfMatrixBirth[id%Nf_][id/Nf_]);
					scalar dec1 = weight*volOwnNode[k]*dfOwnNode[l];
					scalar dec2 = weight*volMatrixDeath[k][id/Nv_]*dfMatrixDeath[id%Nf_][id/Nf_];
					collisionVF_[k][l] +=inc-dec1-dec2;

					if(l<Nf_)
					{
						collisionA_[k][l] += inc*(area_+DatasList[id].area_)-dec1*area_-dec2*DatasList[id].area_;
					}
				}
			}

			//update Matrix
			for(label id=0; id<nNode_; ++id)
			{
				dfMatrixBirth[id%Nf_][id/Nf_] *= dfMatrixBirthForUpdate[id%Nf_][id/Nf_];
				dfMatrixDeath[id%Nf_][id/Nf_] *= DatasList[id].Df_;
			}
		}
	}
	//multiplying constant
	for(label icol=0; icol<nMomentf_; ++icol)
	{
		if(icol ==0)
		{
			for(label irow=0; irow<nMomentv_; ++irow)
			{
				collisionVF_[irow][icol] = collisionVF_[irow][icol]*constantOfCollision;
			}
		}else
		{
			for(label irow=0; irow<Nv_; ++irow)
			{
				collisionVF_[irow][icol] = collisionVF_[irow][icol]*constantOfCollision;
			}
		}
	}

	for(label icol=0; icol<Nv_; ++icol)
	{
		for(label irow=0; irow<Nv_; ++irow)
		{
			collisionA_[irow][icol] = collisionA_[irow][icol]*constantOfCollision;
		}
	}
}


scalar Foam::particleDatas::tauCMin(const List<particleDatas>& DatasList)  //OK
{
	scalar tauCMin(tauC_[0]);
	//compute smallest tauC
	for(label iOwn=0; iOwn<nNode_; ++iOwn)
	{
		for(label iOther=0; iOther<nNode_; ++iOther)
		{
			tauCMin = Foam::min(tauCMin, DatasList[iOwn].tauC_[iOther]);
		}
	}
	return tauCMin;
}

// Busser et al. (2011)
//+POC model , Johannessen et al. (2000 , 2001)
void Foam::particleDatas::sinteringRestructuring(const scalar& deltaTPLoop)  //OK
{
	//fully coalecsed surface area & radius
	am_ = cbrt(M_PI*36.*volume_*volume_);
	scalar rs = 3*volume_/area_;

	//scalar delta =1;

	//characteristic sintering time
	//Eggersdorfer, 2011 + Kirchhof, 2012
	tauS_ = (1-0.367879)*rs*rs*exp(44000/T_)*pow(Np_/2, 0.5-Df_/6)/(1.5874*0.00035);

	//update particle parameter
	if(tauS_ < deltaTPLoop)
	{
		//become sphere
		area_ = am_;
//		k_=1;
		Df_=3;
	}
	else
	{
		//area variation
		scalar sinterd_area = area_-(area_-am_)*deltaTPLoop/tauS_;
		area_ = max(am_,sinterd_area);

		//fractal dimension variation
		//Eggersdorfer, 2011
		scalar c1 = 2.5-0.64*Df_;
		scalar c2 = 1.1*Df_-0.0044*Np_-2.77;
		scalar c3 = 0.015*Np_-6.3*Df_+14.16;

		Df_ = CONST_.Df0 - c1*exp(-sqr((deltaTPLoop/tauS_+c2)/c3));

		//lacunarity variation
		//Lapuerta, 2010 & 2014
		/*とりあえずなし
		scalar Delta = 1 + (0.15*deltaTPLoop)/tauS_;
		scalar c4 = 0.7967*Df_*Df_ -2.38985*Df_ +1.59325;
		scalar c5 = -0.934*Df_*Df_ + 2.86985*Df_ - 0.20385;

		//number of structure
		scalar Ns = (k_ - c4)/c5;

		//compute new k
		k_ = (0.8660254*Ns*Delta - 1.8*Ns*Delta*Delta + 0.7653936*Delta*Delta*Delta)*Df_*Df_
				+ (-1.4433757*Ns*Delta + 7.2*Ns*Delta*Delta -2.2961808*Delta*Delta*Delta)*Df_
					+ (5.1961524*Ns*Delta - 5.4*Ns*Delta*Delta + 1.5307872*Delta*Delta*Delta);
					*/
	}
}



/*　とりあえず抜かしてやってみる
void Foam::particleDatas::surfaceGrowth(const scalar& deltaTPLoop)
{
}
*/


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //
/*
void Foam::particleDatas::operator=(const particleDatas& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::particleDatas::operator=(const Foam::particleDatas&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}
*/
// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //
//weightf_ = weightf[nodeID_][cid_];

// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
