/*
 * particleCellLoop.H
 *
 *  Created on: 2016/05/24
 *      Author: root
 */

nIteration = nIteration*0;

Info << "starting loop" << endl;

forAll(mesh.cells(), cid)
{
	//check message
	if(cid == 1500)
	{
		Info << "Loop is running" << endl;
	}

    // set node data from field(PtrList)
    for(label idatasList=0; idatasList<nNode; ++idatasList)
    {
        datasList[idatasList].setcid(cid, T, mu, lambdaG, deltaDb, rho);
        datasList[idatasList].getWeightsAbscissas(weightv,weightf,abscissav,abscissaf,abscissaa);
    }

    // compute moments
	#include "reproduceMomentsCell.H"
	// initialize source time loop
    deltaTPLeft = runTime.deltaT(); // TODO 宣言 済
    sourceTimeLoop = true; // TODO 宣言 済

    // source term loop between source time loop is true
    while(sourceTimeLoop)
    {
    	// update iteration
        ++nIteration[cid];

        // compute properties
        for(label idatasList=0; idatasList<nNode; ++idatasList)
        {
            datasList[idatasList].set();
        }
        Info << "OK" << endl;
        for(label idatasList=0; idatasList<nNode; ++idatasList)
        {
        	// Calculate collision term against other particels
            datasList[idatasList].updateBetaC(datasList);
        }
        Info << "OK" << endl;

        // compute source time step
        deltaTPLoop.value() = datasList[0].tauCMin(datasList) / nNode;  // 何かで割った方がいい？

        if (deltaTPLoop < deltaTPLeft)
        {
            deltaTPLeft -= deltaTPLoop;
        }
        else
        {
            deltaTPLoop = deltaTPLeft;
            deltaTPLeft.value() = 0;
            // final step decision
            sourceTimeLoop = false;
        }

        // agglomeration
        for(label idatasList=0; idatasList<nNode; ++idatasList)
        {
        	for(label irow=0; irow<nMomentv; ++irow)
        	{
        		momentSetVF[irow][0] += deltaTPLoop.value()*datasList[idatasList].collisionVF_[irow][0];
        	}

        	for(label icol=1; icol<nMomentf; ++icol)
        	{
        		for(label irow=1; irow<Nv; ++irow)
        		{
        			momentSetVF[irow][icol] += deltaTPLoop.value()*datasList[idatasList].collisionVF_[irow][icol];
        		}
        	}

        	for(label icol=0; icol<Nf; ++icol)
        	{
        		for(label irow=0; irow<Nv; ++irow)
        		{
        			momentSetA[irow][icol] += deltaTPLoop.value()*datasList[idatasList].collisionA_[irow][icol];
        		}
        	}
        }
        Info << "OK" << endl;
        // Nucleation
		#include "nucleation.H"

        // Inversion
		#include "MomentInversion.H"
        Info << "OK" << endl;
        // sintering & restructuring
        for(label idatasList=0; idatasList<Nv; ++idatasList)
        {
            datasList[idatasList].set();
            datasList[idatasList].sinteringRestructuring(deltaTPLoop.value());
        }

        // 再生
		#include "reproduceMomentsCell.H"
    }
    //Exit source time loop HERE!!

    // compute properties
    for(label idatasList=0; idatasList<Nv; ++idatasList)
    {
        datasList[idatasList].set();
    }

   ///////OUTPUT DATA////////////////////////////////////////////////////////////////////////
   //node & weight
    //about volume node and weight
    for(label id=0; id<Nv; ++id)
    {
    	label volumeID = Nf*id;
        abscissav[id][cid] = datasList[volumeID].volume_;
        weightv[id][cid] = datasList[volumeID].weightv_;
    }

    for(label idatasList=0; idatasList<nNode; ++idatasList)
    {
        abscissaf[idatasList][cid] = datasList[idatasList].Df_;
        abscissaa[idatasList][cid] = datasList[idatasList].area_;
        weightf[idatasList][cid] = datasList[idatasList].weightf_;
//        weighta[idatasList][cid] = datasList[idatasList].weighta_;
        Dp[idatasList][cid] = datasList[idatasList].Dp_;
//        Dc[idatasList][cid] = datasList[idatasList].Dc_;
        D[idatasList][cid] = datasList[idatasList].D_;
        tauS[idatasList][cid] = datasList[idatasList].tauS_;
        //Lac[idatasList][cid] = datasList[idatasList].k_;
    }

    //moment
    //MOMv
    for(label iMOM=0; iMOM<nMomentv; ++iMOM)
    {
    	MOMv[iMOM][cid] = momentSetVF[iMOM][0];
    }
    //MOMf
    for(label icol=1; icol<nMomentf; ++icol)
    {
    	for(label irow=0; irow<Nf; ++irow)
    	{
    		label iMOM = Nv*(icol-1)+irow;
    		MOMf[iMOM][cid] = momentSetVF[irow][icol];
    	}
    }
    //MOMa
    for(label icol=0; icol<Nf; ++icol)
    {
        for(label irow=0; irow<Nv; ++irow)
        {
        	label iMOM = Nv*icol+irow;
        	MOMf[iMOM][cid] = momentSetA[irow][icol];
        }
     }
    /////OUTPUT END & NEXT CELL///////////////////////////////////////////////////////////////
}
