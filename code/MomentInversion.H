{

//////////////// VOLUME NODES INVERSION //////////////////////////////////////////////////
//output : wi , vi

//Matrix for PD algorithm
label Nv_real = Nv;
scalarSquareMatrix Pv(2*Nv+1, 2*Nv+1, 0);

//set Pv matrix
if (momentSetVF[0][0] < 1.)
{
	//If particle number is less than 1 => There in no Node!!
    Nv_real = 0;
}
else
{
    Pv[0][0] = 1.;

    for(label j=0; j<2*Nv; j+=2)
    {
        Pv[j][1] = momentSetVF[j][0];
    }

    for(label j=1; j<2*Nv; j+=2)
    {
        Pv[j][1] = -momentSetVF[j][0];
    }

    for (label j=2; j<=2*Nv; ++j)
    {
    	//construct recursion
        for (label i=0; i<=2*Nv-j; ++i)
        {
            Pv[i][j] = Pv[0][j-1]*Pv[i+1][j-2]-Pv[0][j-2]*Pv[i+1][j-1];
        }

        if(Pv[0][j] < 1.e-300 || fabs(Pv[0][j]) < fabs(Pv[0][j-1] * Pv[0][j-2])/1.e5 )
        {
            Nv_real = j/2;
            break;
        }
    }
}

//compute jacobian and eigenvalue
if(Nv_real == 0)
{
    for(label i=0; i<nNode; ++i)
    {
    	//no particle => set default value
        datasList[i].volume_ = v_mol;
        datasList[i].Df_ = Df0;
        datasList[i].area_ = a_mol;
        datasList[i].weightv_ = 0.;
        datasList[i].weightf_ = 0.;
    }
}
else
{
	if(Nv_real == 1)
	{
		//nodeID =1,2,3
		for(label i=0; i<Nf; ++i)
		{
			datasList[i].volume_ = momentSetVF[1][0]/stabilise(momentSetVF[0][0], SMALL);
			datasList[i].weightv_ = momentSetVF[0][0];
		}
	}
	else
	{
		//get weight & abscissa of volume
		//construct Jacobian
		scalarList zetav(2*Nv_real, 0.);
		for(label i=1; i<2*Nv_real; ++i)
		{
			if(Pv[0][i]*Pv[0][i-1] > 0.0)
			{
				zetav[i] = Pv[0][i+1]/stabilise(Pv[0][i]*Pv[0][i-1], SMALL);
			}
			else
			{
				zetav[i] = 0.;
			}
		}

		scalarList Aalphav(Nv_real, 0.);
		scalarList Balphav(Nv_real-1, 0.);
		//compute a_alpha
		for(label i=0; i<Nv_real; ++i)
		{
			Aalphav[i] = zetav[2*i]+zetav[2*i+1];
		}
		//compute a_alpha
		for(label i=0; i<Nv_real-1; ++i)
		{
			scalar Balphav_sub = zetav[2*i+1]*zetav[2*i+2];
			Balphav[i] = -1*Foam::sqrt(Balphav_sub);
		}

		//construct Jacobian
		scalarSquareMatrix Jacobianv(Nv_real, Nv_real, 0.);
		for(label i=0; i<Nv_real; ++i)
		{
			Jacobianv[i][i] = Aalphav[i];
		}
		for(label i=0; i<Nv_real-1; ++i)
		{
			Jacobianv[i][i+1] = Balphav[i];
			Jacobianv[i+1][i] = Balphav[i];
		}

		//compute eigen value
		scalarList abscissavValue = eigenValues(Jacobianv);
		//sort eigenvalue as ascending order
		Foam::sort(abscissavValue);
		//delete volume less than v_mol(monomar size)
		/*
		for(label i=0; i<Nv_real; ++i)
		{
			if(abscissavValue[i] < v_mol-1.e-30)
			{
				for(label iValue=0; iValue<Nv_real-1; ++iValue)
		        {
					abscissavValue[i] = abscissavValue[i+1];
		        }
		            --Nv_real;
		            Info << "delete small volume node " << cid << endl;
		     }
		 }
		*/
		//set inversioned data
		for(label i=0; i<Nv_real; ++i)
		{
			scalarList eigenV = eigenVector(Jacobianv, abscissavValue[i]);
			scalar Weight = momentSetVF[0][0]*sqr(eigenV[0]);

			for(label j =0; j<Nf; ++j)
			{
				label idatasList = Nf*i+j;
				datasList[idatasList].volume_ = abscissavValue[i];
				datasList[idatasList].weightv_ = Weight;
			}
		}
	/*
	for(label i=Nv_real; i<Nv; ++i)
	{
		datasList[i].volume_ = datasList[i-1].volume_;
	    datasList[i].weightv_ = 0.;
	}
	*/
//////////////// FRACTAL DIMENSION NODES INVERSION ////////////////////////////////////////////
//output : wij , vij
		//sub Moment for solving fractal dimension
		scalarRectangularMatrix subMomentDf(Nv_real, nMomentf, 0);
		scalarRectangularMatrix volMatrix(Nv_real, Nv_real, 1);
		scalarRectangularMatrix volWeightMatrix(Nv_real, Nv_real, 0);

		//set value
		for(label icol=0; icol<Nv_real; ++icol)
		{
			for(label irow=1; irow<Nv_real; ++irow)
			{
				volMatrix[irow][icol]=volMatrix[irow-1][icol]* datasList[Nf*icol].volume_;
			}
		}
		for(label i=0; i<Nv_real; ++i)
		{
			volWeightMatrix[i][i] = datasList[Nf*i].weightv_;
		}

		//matrix inverse
		scalarRectangularMatrix inverseVolMatrix = SVDinv(volMatrix,0);
		scalarRectangularMatrix inverseVolWeightMatrix = SVDinv(volWeightMatrix,0); //ほんとは不要
		scalarRectangularMatrix inversionMatrix = inverseVolWeightMatrix*inverseVolMatrix;

		//compute subMoments
		for(label icol=0; icol<nMomentf; ++icol)
		{
			scalarList subVector(Nv_real, 0.);
			scalarList subMomentVector(Nv_real, 0.);

			//set subvector data
			for(label i=0; i<Nv_real; ++i)
			{
				subVector[i] = momentSetVF[i][icol];
			}

			//compute subMomentvector
			for(label i=0; i<Nv_real; ++i)
			{
				for(label j=0; j<Nv_real; ++j)
				{
					subMomentVector[i] += inversionMatrix[i][j]*subVector[j];
				}
			}

			//set to subMomentDf
			for(label irow=0; irow<Nv_real; ++irow)
			{
				subMomentDf[irow][icol] = subMomentVector[irow];
			}
		}

		//PD algorithm for solving nodes aboud Df each volume
		for(label volumeID=0; volumeID<Nv_real; ++volumeID)
		{
			//construct matrix P
			scalarSquareMatrix PDf(2*Nf+1, 2*Nf+1, 0);
			//set data to P
			PDf[0][0] = 1.;
			for(label j=0; j<2*Nf; j+=2)
			{
				PDf[j][1] = subMomentDf[volumeID][j];
			}
			for(label j=1; j<2*Nf; j+=2)
			{
				PDf[j][1] = -subMomentDf[volumeID][j];
			}

			for (label j=2; j<=2*Nf; ++j)
			{
				//construct recursion
				for (label i=0; i<=2*Nf-j; ++i)
				{
					PDf[i][j] = PDf[0][j-1]*PDf[i+1][j-2]-PDf[0][j-2]*PDf[i+1][j-1];
				}
			}

			//compute jacobian
			scalarList zetaDf(nMomentf, 0.);
			for(label i=1; i<nMomentf; ++i)
			{
				if(PDf[0][i]*PDf[0][i-1] > 0.)
				{
					zetaDf[i] = PDf[0][i+1]/stabilise(PDf[0][i]*PDf[0][i-1], SMALL);
				}
				else
				{
					zetaDf[i] = 0.;
				}
			}
			scalarList AalphaDf(Nf, 0.);
			scalarList BalphaDf(Nf, 0.);
			//compute a_alpha
			for(label i=0; i<Nf; ++i)
			{
				AalphaDf[i] = zetaDf[2*i]+zetaDf[2*i+1];
			}
			//compute a_alpha
			for(label i=0; i<Nf-1; ++i)
			{
				scalar BalphaDf_sub = zetaDf[2*i+1]*zetaDf[2*i+2];
				BalphaDf[i] = - Foam::sqrt(BalphaDf_sub);
			}

			//construct Jacobian
			scalarSquareMatrix JacobianDf(Nf, Nf, 0.);
			for(label i=0; i<Nf; ++i)
			{
				JacobianDf[i][i] = AalphaDf[i];
			}
			for(label i=0; i<Nf-1; ++i)
			{
				JacobianDf[i][i+1] = BalphaDf[i];
				JacobianDf[i+1][i] = BalphaDf[i];
			}

			//compute eigen value
			scalarList abscissaDfValue = eigenValues(JacobianDf);
			//sort eigenvalue
			Foam::sort(abscissaDfValue);

			/*
			 *ALGORITHM TO CORRECT VALUE  but is this OK??
			 */
			for(label i=0; i<Nf; ++i)
			{
				if(abscissaDfValue[i] < 1.)
				{
					abscissaDfValue[i]=1.;
				}
				else if(abscissaDfValue[i] > 3.)
				{
					abscissaDfValue[i]=3.;
				}
			}

			//set inversioned data
			for(label i=0; i<Nf; ++i)
			{
				datasList[Nf*volumeID+i].Df_ = abscissaDfValue[i];
				datasList[Nf*volumeID+i].weightf_ = sqr(eigenVector(JacobianDf, abscissaDfValue[i])[0]);
			}
		}

//////////////// AREA NODES INVERSION ///////////////////////////////////////////////////////
//output : aij
		scalarRectangularMatrix subMomentA(Nv_real, Nf, 0);
		scalarRectangularMatrix subsubMomentA(Nv_real, Nf, 0);

		//set data to subMomentA
		for(label icol=0; icol<Nf; ++icol)
		{
			for(label irow=0; irow<Nv_real; ++irow)
			{
				for(label i=0; i<Nv_real; ++i)
				{
					subMomentA[irow][icol] += inversionMatrix[irow][i]*momentSetA[i][icol];
				}
			}
		}

		//compute subsubMomentA of each volume node
		for(label volumeID=0; volumeID<Nv_real; ++volumeID)
		{
			scalarRectangularMatrix dfMatrix(Nf, Nf, 1);
			scalarRectangularMatrix dfWeightMatrix(Nf, Nf, 0);

			//set value
			for(label icol=0; icol<Nf; ++icol)
			{
				scalar df_onV1 = datasList[Nf*volumeID+icol].Df_;
				for(label irow=1; irow<Nf; ++irow)
				{
					dfMatrix[irow][icol]=dfMatrix[irow-1][icol]* df_onV1;
				}
			}
			for(label icol=0; icol<Nf; ++icol)
			{
				dfWeightMatrix[icol][icol] = datasList[Nf*volumeID+icol].weightf_;
			}

			//compute inversion matrix
			scalarRectangularMatrix inverseDfMatrix = SVDinv(dfMatrix,0);
			scalarRectangularMatrix inverseDfWeightMatrix = SVDinv(dfWeightMatrix,0); //ほんとは不要
			scalarRectangularMatrix inversionDfMatrix = inverseDfWeightMatrix * inverseDfMatrix;

			//compute subsubMomentA
			for(label icol=0; icol<Nf; ++icol)
			{
				for(label i=0; i<Nf; ++i)
				{
					subsubMomentA[volumeID][icol] += inversionDfMatrix[icol][i]*subMomentA[icol][i];
				}
			}
		}

		//area node inversion
		for(label i=0; i<Nf; ++i)
		{
			for(label j=0; j<Nv_real; ++j)
			{
				label id = Nf*i+j;
				datasList[id].area_ = subsubMomentA[i][j];
			}
		}
	}

}
}
