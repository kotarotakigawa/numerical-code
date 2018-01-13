/*
 * nucleation.H
 *
 *  Created on: 2017/12/06
 *      Author: root
 */

#ifndef NUCLEATION_H_
#define NUCLEATION_H_

	 //Nucleation//////////////////////////////////////////////////////////////////

	scalar sutRatioI = sutRatio[cid];

	if(sutRatioI>1)
	{
		//Info << "starting nucleation cid" << cid << endl;
		scalar nucRateI = nucRate[cid];
		scalar nucDiamI = nucDiam[cid];
		scalar SIO2GI = SIO2G[cid];
		scalar rhoI = rho[cid];

        scalar nucVolume = M_PI*nucDiamI*nucDiamI*nucDiamI/6;
        scalar nucSurface = M_PI*nucDiamI*nucDiamI;
        scalar nucNumber = min(nucRateI*deltaTPLoop.value(), SIO2GI*rhoI/(M_p*nucVolume));
        nucNumber = max(0, nucNumber);


        //prepare moment matrix
        scalarRectangularMatrix nucleationVF(nMomentv, nMomentf, 0.);
        scalarRectangularMatrix nucleationA(Nv, Nf, 0.);
        nucleationVF[0][0] = nucNumber;
        nucleationA[0][0] = nucNumber*nucSurface;

        // calculate variation of SIO2 concentration by nucleation
        for(label irow=1; irow<nMomentv; ++irow)
        {
        	nucleationVF[irow][0] = nucVolume*nucleationVF[irow-1][0];
        }
        for(label irow=0; irow<Nv; ++irow)
        {
        	for(label icol=1; icol<nMomentf; ++icol)
        	{
        		nucleationVF[irow][icol] = nucleationVF[irow][icol-1]*Df0;
        	}
        }
        for(label irow=1; irow<Nv; ++irow)
        {
        	nucleationA[irow][0] = nucleationA[irow-1][0]*nucVolume;
        }
        for(label irow=0; irow<Nv; ++irow)
        {
        	for(label icol=1; icol<Nf; ++icol)
        	{
        		nucleationA[irow][icol] = nucleationA[irow][icol-1]*Df0;
        	}
        }

        //OUTPUT to momentSet
        //apply moment variation

        	for(label irow=0; irow<nMomentv; ++irow)
        	{
        		for(label icol=0; icol<nMomentf; ++icol)
        		{
        			momentSetVF[irow][icol] += nucleationVF[irow][icol];
        		}
        	}

        	for(label icol=0; icol<Nf; ++icol)
        	{
        		for(label irow=0; irow<Nv; ++irow)
        		{
        			momentSetA[irow][icol] += nucleationA[irow][icol];
        		}
        	}
	}

#endif /* NUCLEATION_H_ */
