/*
 * reproduceMomentsCell.H
 *
 *  Created on: 2016/05/23
 *      Author: root
 */

//initialize
for(label i=0; i<nMomentv; ++i)
{
	for(label j=0; j<nMomentf; ++j)
	{
		momentSetVF[i][j] = 0.;
	}
}

for(label i=0; i<Nv; ++i)
{
	for(label j=0; j<Nf; ++j)
	{
		momentSetA[i][j] = 0.;
	}
}

//prepare matrix & vectors
//volume matrix
//(row, column)
scalarRectangularMatrix Vm(nMomentv, Nv , 1.);
//volume weight vector
scalarList Vwv(Nv, 0.);
//set data
for(label icol = 0; icol < Nv; ++icol)
{
	label iData = icol/Nv;
	Vwv[icol] = datasList[iData].weightv_;

	for(label irow = 1; irow < nMomentv; ++irow)
	{
		Vm[irow][icol] = Vm[irow - 1][icol]*datasList[iData].volume_;
	}
}
//volume weighted square matrix for fractal dimension
scalarSquareMatrix Vwm(Nv, Nv , 0.);
for(label irow = 0; irow < Nv; ++irow)
{
	for(label icol = 0; icol < Nv; ++icol)
	{
		Vwm[irow][icol] = Vm[irow][icol]*Vwv[icol];
	}
}

//fractal weight matrix
scalarSquareMatrix Fwm(Nv, Nf, 0);
//fractal dimension matrix
scalarSquareMatrix Fdm(Nv, Nf, 0);
//weighted fractal matrix for Loop : wijDfij^l
scalarSquareMatrix FwmLoop(Nv, Nf, 0);
//set data
for(label irow =0; irow < Nv; ++irow){
	for(label icol =0; icol < Nf; ++icol)
	{
		Fwm[irow][icol] = datasList[Nf*irow+icol].weightf_;
		Fdm[irow][icol] = datasList[Nf*irow+icol].Df_;
		FwmLoop[irow][icol] = Fwm[irow][icol]*Fdm[irow][icol];
	}
}

//area matrix
scalarSquareMatrix Am(Nv, Nf, 0);
//weighted fractal matrix for Loop : wijaij^l
scalarSquareMatrix AwmLoop(Nv, Nf, 0);
//setdata
for(label irow =0; irow < Nv; ++irow){
	for(label icol =0; icol < Nf; ++icol)
	{
		Am[irow][icol] = datasList[Nf*irow+icol].area_;
		AwmLoop[irow][icol] = Am[irow][icol]*Fwm[irow][icol];
	}
}


//produce moment matrix (mk00 : volume only)
for(label irow =0; irow < nMomentv; ++irow){
	for(label i =0; i < Nv; ++i)
	{
		momentSetVF[irow][0] += Vm[irow][i]*Vwv[i];
	}
}


//produce moment matrix (mkl0 : volume & fractal dimension)
for(label icol =1; icol < nMomentf; ++icol)
{
	for(label irow =0; irow < Nv; ++irow)
	{
		scalarList FwvLoop(Nv, 0.);
		//calculate FwvLoop
		for(label i =0; i < Nv; ++i)
		{
			for(label j=0; j<Nf; ++j)
			{
				FwvLoop[i] += FwmLoop[i][j];
			}
		}
		//calculate momentSet
		for(label i =0; i < Nv; ++i)
		{
			momentSetVF[irow][icol] += Vwm[irow][i]*FwvLoop[i];
		}
	}

	//update FwmLoop
	for(label i =0; i < Nv; ++i)
	{
		for(label j=0; j<Nf; ++j)
		{
			FwmLoop[i][j] *= Fdm[i][j];
		}
	}
}


//produce moment matrix (mkl1 : volume & fractal dimension & surface area)
for(label icol =1; icol < Nv; ++icol)
{
	for(label irow =0; irow < Nf; ++irow)
	{
		scalarList AwvLoop(Nv, 0.);
		//calculate AwLoop
		for(label i =0; i < Nv; ++i)
		{
			for(label j=0; j<Nf; ++j)
			{
				AwvLoop[i] += AwmLoop[i][j];
			}
		}
		//calculate momentSet
		for(label i =0; i < Nv; ++i)
		{
			momentSetA[irow][icol] += Vwm[irow][i]*AwvLoop[i];
		}
	}

	//update FwmLoop
	for(label i =0; i < Nv; ++i)
	{
		for(label j=0; j<Nf; ++j)
		{
			AwmLoop[i][j] *= Fdm[i][j];
		}
	}
}

