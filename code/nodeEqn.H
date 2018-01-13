//nodeEqn.H

volScalarField weightedVolLoop("weightedVolLoop", weightv[0]);
//initialize weightedVolLoop

forAll(mesh.cells(), cid)
{
	weightedVolLoop[cid] = 0.;
}

for(label iNode=0; iNode<nNode; ++iNode)
{
	label volumeID = iNode/Nf;

	////WiWij
	volScalarField weightI("weightI", weightv[volumeID]*weightf[iNode]);
	surfaceScalarField weightISurf("weightISurf", fvc::interpolate(weightI));

	//weighted scalar(center of cell)
	volScalarField weightedVol("weightedVol", weightI*abscissav[volumeID]);
	volScalarField weightedDf("weightedDf", weightI*abscissaf[iNode]);
	volScalarField weightedShape("weightedShape", weightI*shape[iNode]);

	//weighted scalar(surface of cell)
	surfaceScalarField weightedVolSurf("weightedVolSurf", fvc::interpolate(weightedVol));
	surfaceScalarField weightedDfSurf("weightedDfSurf", fvc::interpolate(weightedDf));
	surfaceScalarField weightedShapeSurf("weightedShapeSurf", fvc::interpolate(weightedShape));

	//solve eqs
	solve
	(
			fvm::ddt(weightI)
				+fvc::div(weightISurf*UphiP)
				-fvc::laplacian(D[iNode],weightI)
	);
	solve
	(
			fvm::ddt(weightedVol)
				+fvc::div(UphiP*weightedVolSurf)
				-fvc::laplacian(D[iNode],weightedVol)
	);
	solve
	(
			fvm::ddt(weightedDf)
   				+fvc::div(UphiP*weightedDfSurf)
   				-fvc::laplacian(D[iNode],weightedDf)
	);
	solve
	(
			fvm::ddt(weightedShape)
        		+fvc::div(UphiP*weightedShapeSurf)
        		-fvc::laplacian(D[iNode],weightedShape)
	);


	//provisional insert & correction after
	forAll(mesh.cells(), cid)
	{
		//about volume
		weightedVolLoop[cid] += weightedVol[cid];
		weightf[iNode][cid] = weightI[cid];
		//about Df & area
		abscissaf[iNode][cid] = weightedDf[cid];
		shape[iNode][cid] = weightedShape[cid];
	}

	//Wi calculation(when iNode=2,5,8) & correction
	if(iNode%Nf == Nf-1)
	{
		forAll(mesh.cells(), cid)
		{
			//calculate weightSum
			scalar weightSum = 0.;
			for(label i=0; i<Nf; ++i)
			{
				weightSum += weightf[iNode-i][cid];
			}
			//set value Wi
			weightv[volumeID][cid] = weightSum;

			//set value Wij
			for(label i=0; i<Nf; ++i)
			{
				if(weightSum <= 0.)
				{
					weightf[iNode-i][cid] = 0.;
				}
				else
				{
					weightf[iNode-i][cid] /= weightSum;
				}
			}
		}

		//correction
		//about volume abscissa
		forAll(mesh.cells(), cid)
		{
			if(weightv[volumeID][cid] < 1.)
			{
				abscissav[volumeID][cid] = v_mol;
			}
			else
			{
				abscissav[volumeID][cid] = weightedVolLoop[cid]/weightv[volumeID][cid];
			}
		}

		//about other scalar
		for(label i=0; i<Nf; ++i)
		{
			label id = iNode - i;

			forAll(mesh.cells(), cid)
			{
				if(weightv[volumeID][cid] < 1.)
				{
					weightv[volumeID][cid] = 0.;
					weightf[id][cid] = 0.;
					shape[id][cid] = 1.;
					abscissaa[id][cid] = a_mol;
					abscissaf[id][cid] = Df0;
				}
				else
				{
					shape[id][cid] = max(0., shape[id][cid]/(weightv[volumeID][cid]*weightf[id][cid]));
					shape[id][cid] = min(1., shape[id][cid]);
					abscissaa[id][cid] = max(a_mol, a_mol*Foam::pow(abscissav[volumeID][cid]/v_mol, 1-shape[id][cid]/3));
					abscissaf[id][cid] = max(1., abscissaf[id][cid]/(weightv[volumeID][cid]*weightf[id][cid]));
				}
			}
		}
		//reset weightedVolLoop
		forAll(mesh.cells(), cid)
		{
			weightedVolLoop[cid] = 0.;
		}
	}
}
//TODO 吐き出しの調整
