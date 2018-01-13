/*
 * createQMOM.H
 *
 *  Created on: 2015/05/15
 *      Author: Matsushita H
 */

#ifndef QMOMFOAM_CREATEQMOM_H_
#define QMOMFOAM_CREATEQMOM_H_

	Info << "Creating constants" << endl;
	//create constants
    const label Nv(particleDatas::CONST_.Nv); // vノード数　可変
    const label Nf(particleDatas::CONST_.Nf); // Dfノード数　可変
    const scalar Df0(particleDatas::CONST_.Df0);
    const scalar kf(particleDatas::CONST_.kf);

    //number of moments
	const label nMomentv(particleDatas::CONST_.nMomentv); // Radau-CQMOM用
	const label nMomentf(particleDatas::CONST_.nMomentf);
	const label nMoment(particleDatas::CONST_.nMoment); // Radau-CQMOM用
	const label nNode(particleDatas::CONST_.nNode);

	// node information for cell Loop
	/*
	 * EXAMPLE
	 * i=0 v1 Df11 // i=3 v2 Df21 // i=6 v3 Df31
	 * i=1 v1 Df12 // i=4 v2 Df22 // i=7 v3 Df32
	 * i=2 v1 Df13 // i=5 v2 Df23 // i=8 v3 Df33
	 */
	List<particleDatas> datasList(nNode);

	//Prepare momentsets for cell Loop
    scalarRectangularMatrix momentSetVF(nMomentv, nMomentf, 0);  //volume and fractal dimension
    scalarRectangularMatrix momentSetA(Nv, Nf, 0);  //surface area

	//Boltzmann constant
    const dimensionedScalar kb("kb", dimensionSet(1,2,-2,-1,0,0,0), 1.38e-23);
    //gas constant
    const dimensionedScalar Rg("Rg", dimensionSet(1,2,-2,-1,-1,0,0), 8.31);
    //???????
    const dimensionedScalar Ea("Ea", dimensionSet(1,2,-2,0,-1,0,0), 4.1e5);
    //SiO2 monomar's volume
    const scalar v_mol(particleDatas::CONST_.v_mol);
    //SiO2 monomar's surface area
    const scalar a_mol(particleDatas::CONST_.a_mol);
    //SiO2 monomar's diameter
    const scalar d_mol(particleDatas::CONST_.d_mol);
    //diameter of air molecular
    const dimensionedScalar da("da", dimensionSet(0,1,0,0,0,0,0), 3.70e-10);
    //avogadro number
    const dimensionedScalar Nav("Nav", dimensionSet(0,0,0,0,-1,0,0), 6.022e23);
    //????????
    const dimensionedScalar deltaD0("deltaD0", dimensionSet(0,3,-1,0,0,0,0), 2.36e-12);

    const scalar Df_asymp(particleDatas::CONST_.Df_asymp);
    const volScalarField deltaDb = deltaD0*exp(-Ea/(Rg*T));
    const volScalarField lambdaG = kb*T/Foam::sqrt(2.)/M_PI/p/pow(da,2);

    //set time controll
    // Left time for cell Loop (time step written in controllDict)
    dimensionedScalar deltaTPLeft(runTime.deltaT());
    // Time step for cell Loop
    dimensionedScalar deltaTPLoop(deltaTPLeft);
    // decision boolean for cell Loop do /or not?
    bool sourceTimeLoop;


    Info << "Creating variables" << endl;

    //Prepare MOM framework only
    PtrList<volScalarField> MOMv(nMomentv);
    //Prepare MOMk00 (weighted by volume only)
    for(label id = 0; id < nMomentv; ++id)
    {
        std::stringstream   ssMOMv;

        ssMOMv << "MOMv" << id << 0  << 0;

        word wMOMv = ssMOMv.str();

        Info << "Reading field  " << wMOMv << endl;
        MOMv.set
        (
                id,
                new volScalarField
                (
                        IOobject
                        (
                                wMOMv,
                                runTime.timeName(),
                                mesh,
                                IOobject::READ_IF_PRESENT,
                                IOobject::AUTO_WRITE
                        ),
                        mesh,
                        dimensionedScalar(wMOMv, dimensionSet(0, 3*id-3, 0, 0, 0, 0, 0), 0.0)
                )
        );
    }

    //Prepare MOMkl0 (weighted by volume and fractal dimension)
    //excepting MOMi00
    PtrList<volScalarField> MOMf(nMomentf*Nv-Nv);
    for(label id = 0; id < nMomentf*Nv-Nv; ++id )
   	{
    	std::stringstream   ssMOMf;
        ssMOMf << "MOMf" << id%Nv << id/Nv+1 << 0;

        word wMOMf = ssMOMf.str();

        Info << "Reading field " << wMOMf << endl;
        MOMf.set
        (
        	id,
        	new volScalarField
        	(
        		IOobject
        		(
        			wMOMf,
        			runTime.timeName(),
        			mesh,
        			IOobject::READ_IF_PRESENT,
        			IOobject::AUTO_WRITE
        		),
        		mesh,
        		dimensionedScalar(wMOMf, dimensionSet(0, 3*(id%Nv)-3, 0, 0, 0, 0, 0), 0.0)
        	)
        );
   	}

    //Prepare MOMkl1 (weighted by volume, fractal dimension and surface area)
    PtrList<volScalarField> MOMa(nNode);
    for(label id = 0; id < nNode; ++id )
    {
    	std::stringstream   ssMOMa;
    	ssMOMa << "MOMa" << id%Nf << id/Nf << 1;
   		word wMOMa = ssMOMa.str();

    	Info << "Reading field " << wMOMa << endl;
   		MOMa.set
   		(
  			id,
  			new volScalarField
    		(
   				IOobject
   				(
   					wMOMa,
   					runTime.timeName(),
   					mesh,
   					IOobject::READ_IF_PRESENT,
   					IOobject::AUTO_WRITE
    			),
    			mesh,
   				dimensionedScalar(wMOMa, dimensionSet(0, 3*(id%Nf)+2-3, 0, 0, 0, 0, 0), 0.0)
   			)
    	);
    }



    //Prepare volume weight
    PtrList<volScalarField> weightv(Nv);

    for(label id = 0; id < Nv; ++id)
    {
        std::stringstream   ssweightv;

        ssweightv << "w" << id+1 ;

        word wweightv = ssweightv.str();

        Info << "Reading field " << wweightv << endl;
        weightv.set
        (
        	id,
            new volScalarField
            (
            	IOobject
                (
                	wweightv,
                	runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                 ),
                 mesh,
                 dimensionedScalar(wweightv, dimensionSet(0, -3, 0, 0, 0, 0, 0), 0.)
             )
        );
    }

    //Prepare volume abcissa
    PtrList<volScalarField> abscissav(Nv);

    for(label id = 0; id < Nv; ++id)
    {
        std::stringstream   ssabscissav;

        ssabscissav << "v" << id+1 ;

        word wabscissav = ssabscissav.str();

        Info << "Reading field " << wabscissav << endl;
        abscissav.set
        (
        	id,
            new volScalarField
            (
            	IOobject
                (
                	wabscissav,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                 ),
                 mesh,
                 dimensionedScalar(wabscissav, dimensionSet(0, 3, 0, 0, 0, 0, 0), v_mol)
             )
        );
    }

    //Prepare fractal dimension weight
    PtrList<volScalarField> weightf(nNode);
    for(label id = 0; id < nNode; ++id)
    {
    	std::stringstream   ssweightf;

    	ssweightf << "wf" << id/Nf+1 << id%Nf+1 ;

    	word wweightf = ssweightf.str();

    	Info << "Reading field " << wweightf << endl;
    	weightf.set
    	(
    		id,
    		new volScalarField
    		(
    			IOobject
    			(
              		wweightf,
              		runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
    			mesh,
    			dimensionedScalar(wweightf, dimless, 0.)
    		)
    	);
    }

    //Prepare fractal dimension abcissa
    PtrList<volScalarField> abscissaf(nNode);

    for(label id = 0; id < nNode; ++id)
   	{
   		std::stringstream   ssabscissaf;

    	ssabscissaf << "Df" << id/Nf+1 << id%Nf+1 ;

   		word wabscissaf = ssabscissaf.str();

   		Info << "Reading field " << wabscissaf << endl;
   		abscissaf.set
   		(
   			id,
            new volScalarField
            (
           		IOobject
           		(
             		wabscissaf,
             		runTime.timeName(),
             		mesh,
             		IOobject::READ_IF_PRESENT,
             		IOobject::AUTO_WRITE
             	),
             	mesh,
             	dimensionedScalar(wabscissaf, dimless, 3.0)  //要修正
             )
   		);
   	}

    //Prepare surface area abcissa
    PtrList<volScalarField> abscissaa(nNode);

    for(label id = 0; id < nNode; ++id)
    {
    	std::stringstream   ssabscissaa;

    	ssabscissaa << "a" << id/Nf+1 << id%Nf+1 ;

   		word wabscissaa = ssabscissaa.str();

   		Info << "Reading field " << wabscissaa << endl;
   		abscissaa.set
   		(
   			id,
    		new volScalarField
    		(
    			IOobject
    			(
    				wabscissaa,
    				runTime.timeName(),
    				mesh,
    				IOobject::READ_IF_PRESENT,
    				IOobject::AUTO_WRITE
    			),
    			mesh,
    			dimensionedScalar(wabscissaa, dimArea, a_mol) //要修正
    		)
    	);
    }

    PtrList<volScalarField> shape(nNode);

    for(label id = 0; id < nNode; ++id)
    {
    	std::stringstream   ssshape;

    	ssshape << "shape" << id/Nf+1 << id%Nf+1;

    	word wshape = ssshape.str();

    	Info << "Reading field " << wshape <<  endl;
    	shape.set
    	(
    		id,
    		new volScalarField
    		(
    			IOobject
    			(
    				wshape,
    				runTime.timeName(),
    				mesh,
    				IOobject::READ_IF_PRESENT,
    				IOobject::AUTO_WRITE
    			),
    			mesh,
    			dimensionedScalar(wshape, dimless, 1.)
    		)
    	);
    }
/*
    volScalarField SIO
    (
    	IOobject
    	(
    		"SIO",minNucDiam
    		runTime.timeName(),
    		mesh,
    		IOobject::MUST_READ,
    		IOobject::NO_WRITE
    	),
    	mesh
    );

    Info<< "Reading field H2O" << endl;
    volScalarField H2O
    (
    	IOobject
    	(
    		"H2O",
    		runTime.timeName(),
    		mesh,
    		IOobject::MUST_READ,
    		IOobject::NO_WRITE
    	),
    	mesh
    );
*/
    // weight fraction
    Info<< "Reading field SIO2G" << endl;
    volScalarField SIO2G
    (
    	IOobject
    	(
    		"SIO2G",
    		runTime.timeName(),
    		mesh,
    		IOobject::MUST_READ,
    		IOobject::AUTO_WRITE
    	),
    	mesh
    );
    SIO2G.max(0.);

    ///////////constants for createField of nucleation/////////////////////
    //prefactor of vapor presure of SIO2
    const dimensionedScalar p0("p0", dimensionSet(1,-1,-2,0,0,0,0), 133.322*3.0e-6);
    //SiO2 weight per volume
    const dimensionedScalar rho_p("rho_p", dimensionSet(1,-3,0,0,0,0,0),particleDatas::CONST_.rho_p);
    //SiO2 weight per mole
    const dimensionedScalar Mp("Mp",dimensionSet(1,0,0,0,-1,0,0), particleDatas::CONST_.M_p);
    const scalar M_p(particleDatas::CONST_.M_p);
    //moleculer volume of SIO2  >>>> this is equal to v_mol
    const dimensionedScalar v0("v0",dimensionSet(0,3,0,0,0,0,0), particleDatas::CONST_.v_mol);
    //moleculer volume of SIO2  >>>> this is equal to a_mol
    const dimensionedScalar a0("a0",dimensionSet(0,2,0,0,0,0,0), particleDatas::CONST_.a_mol);
    //moleculer weight of SIO2
    const dimensionedScalar mp("mp", dimensionSet(1,0,0,0,0,0,0), 1.203e-25);
    //prefactor of vapor pressure of sio2
    const dimensionedScalar preT("preT", dimensionSet(0,0,0,-1,0,0,0), 0.0087);
    //surface tension of SIO2
    const dimensionedScalar sTension("sTension", dimensionSet(1,0,-2,0,0,0,0), 0.307);

    //float
    //const dimensionedScalar Notzero("Notzero", dimensionSet(0,0,0,0,0,0,0), 1e-20);

    Info<< "Calculating field Saturation ratio" << endl;
    volScalarField sutRatio
    (
    	IOobject
    	(
    		"sutRatio",
    		runTime.timeName(),
    		mesh,
    		IOobject::NO_READ,
    		IOobject::AUTO_WRITE
    	),
    	(rho*SIO2G*Nav/Mp)/((p0*exp(preT*T)/(kb*T)))+SMALL
    );
    sutRatio.max(SMALL);

    Info<< "Calculating field Nucleation rate" << endl;
    volScalarField nucRate
    (
    	IOobject
    	(
    		"nucRate",
    		runTime.timeName(),
    		mesh,
    		IOobject::NO_READ,
    		IOobject::AUTO_WRITE
    	),
    	pow(p0*exp(preT*T)/(kb*T), 2)*sutRatio*v0*pow(2*sTension/(M_PI*mp), 1/2)
    		*exp(a0*sTension/(kb*T) - 4*pow(a0*sTension/(kb*T), 3)/(SMALL+27*log(sutRatio)*log(sutRatio)))
    );
    nucRate.max(0.);

    Info<< "Calculating field Nuclei size" << endl;
    volScalarField nucDiam
    (
    	IOobject
    	(
    		"nucDiam",
    		runTime.timeName(),
    		mesh,
    		IOobject::NO_READ,
    		IOobject::AUTO_WRITE
    	),
    	4*sTension*v0/(kb*T*log(sutRatio))
    );
    nucDiam.max(0.);
/*
    Info << "Calculating the smallest nuclei" << endl;
	dimensionedScalar minNucDiam("minNucDiam", dimensionSet(0,1,0,0,0,0,0), 1);
    forALL(mesh.cells(), cid)
    {
    	if(sutRatio[cid]>1)
    	{
    		minNucDiam = min(nucDiam[cid], minNucDiam);
    	}
    }
    const dimensionedScalar minNucVol("minNucVol", dimensionSet(0,3,0,0,0,0,0), M_PI*minNucDiam*minNucDiam*minNucDiam/6);
    const dimensionedScalar minNucArea("minNucArea", dimensionSet(0,2,0,0,0,0,0), M_PI*minNucDiam*minNucDiam);
*/


/*
    const dimensionedScalar AReact("AReact", dimensionSet(0,3,-1,0,-1,0,0), 8.5e4);
    const dimensionedScalar EReact("EReact", dimensionSet(1,2,-2,0,-1,0,0), 5650);
    const dimensionedScalar MwH2O("MwH2O", dimensionSet(1,0,0,0,-1,0,0), 0.018015);
    const dimensionedScalar MwSIO("MwSIO", dimensionSet(1,0,0,0,-1,0,0), 0.044085);

    Info<< "Calculating field SIO2gen" << endl;
    volScalarField SIO2gen
    (
        IOobject
        (
            "SIO2gen",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        AReact*exp(-EReact/Rg/T)*rho*SIO/MwSIO*rho*H2O/MwH2O
    );
*/
/*
    //どうしようか
    ///Moment gneration ==> NEED TO REVISE///
    scalarRectangularMatrix momentGenVF(nMomentv, nMomentf, 0);
    for(label i = 0; i < nMomentv; ++i)
    {
        momentGenVF[i] = pow(V0, i);
    }
    for(label i = nMomentv; i < nMoment ; ++i)
    {
        label k = (i+1)%Nv;
        label l = (i+1)/Nv-1;
        momentGen[i] = pow(V0, k)*pow(A0, l);
    }

    scalarSquareMatrix momentGenA(Nv, Nf, 0);
    for(label i = 0; i < nMomentv; ++i)
    {
    	momentGen[i] = pow(V0, i);
    }
    for(label i = nMomentv; i < nMoment ; ++i)
    {
    	label k = (i+1)%Nv;
    	label l = (i+1)/Nv-1;
    	momentGen[i] = pow(V0, k)*pow(A0, l);
    }
*/
    //////////////////////////////////////////////////

    PtrList<volScalarField> Dp(nNode);

    for(label id = 0; id < nNode; ++id)
   	{
    	std::stringstream   ssDp;

    	ssDp << "Dp" << id/Nf+1 << id%Nf+1 ;

   		word wDp = ssDp.str();

   		Info << "Reading field " << wDp << endl;
   		Dp.set
   		(
   			id,
    		new volScalarField
    		(
 				IOobject
 				(
   					wDp,
   					runTime.timeName(),
    				mesh,
    				IOobject::NO_READ,
    				IOobject::AUTO_WRITE
    			),
    			mesh,
               	dimensionedScalar(wDp, dimLength, d_mol)
    			)
   		);
    }

/*
    PtrList<volScalarField> Vp(nNode);

    for(label id = 0; id < nNode; ++id)
   	{
    	std::stringstream   ssVp;

    	ssVp << "Vp" << id/Nf+1 << id%Nf+1 ;

   		word wVp = ssVp.str();

   		Info << "Reading field " << wVp << endl;
   		Vp.set
   		(
   			id,
    		new volScalarField
    		(
 				IOobject
 				(
   					wVp,
   					runTime.timeName(),
    				mesh,
    				IOobject::NO_READ,
    				IOobject::NO_WRITE
    			),
    			mesh,
               	dimensionedScalar(wVp, dimVolume, v_mol)
    			)
   		);
    }


    PtrList<volScalarField> Ap(nNode);

    for(label id = 0; id < nNode; ++id)
   	{
    	std::stringstream   ssAp;

    	ssAp << "Ap" << id/Nf+1 << id%Nf+1 ;

   		word wAp = ssAp.str();

   		Info << "Reading field " << wAp << endl;
   		Ap.set
   		(
   			id,
    		new volScalarField
    		(
 				IOobject
 				(
   					wAp,
   					runTime.timeName(),
    				mesh,
    				IOobject::NO_READ,
    				IOobject::AUTO_WRITE
    			),
    			mesh,
               	dimensionedScalar(wAp, dimArea, a_mol)
    			)
   		);
    }

    //lacunarity(fractal prefactor)
    PtrList<volScalarField> Lac(nNode);

    for(label id = 0; id < nNode; ++id )
    {
       	std::stringstream   ssLac;
       	ssLac << "Lac" << id/Nf+1 << id%Nf+1 ;

       	word wLac = ssLac.str();

       	Info << "Reading field Lac" << id/Nf+1 << id%Nf+1 <<  endl;
       	Lac.set
       	(
    		id,
    		new volScalarField
    		(
    			IOobject
    			(
    				wLac,
    				runTime.timeName(),
    				mesh,
    				IOobject::NO_READ,
    				IOobject::AUTO_WRITE
    			),
    			mesh,
    			dimensionedScalar(wLac, dimless, 1.0)
    		)
       	);
    }
*/
    PtrList<volScalarField> D(nNode);

    for(label id =0; id < nNode; ++id)
   	{
   		std::stringstream   ssD;

   		ssD << "D" << id/Nf+1 << id%Nf+1 ;

   		word wD = ssD.str();

   		Info << "Reading field " << wD << endl;
    	D.set
   		(
   			id,
    		new volScalarField
    		(
				IOobject
				(
					wD,
					runTime.timeName(),
    				mesh,
    				IOobject::NO_READ,
    				IOobject::AUTO_WRITE
    			),
    			mesh,
    			dimensionedScalar(wD, dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0)
    		)
    	);
   	}

    PtrList<volScalarField> tauS(nNode); //all charactreristic time is needed?

    for(label id = 0; id < nNode; ++id)
   	{
   		std::stringstream   sstauS;

   		sstauS << "tauS" << id/Nf+1 << id%Nf+1 ;

   		word wtauS = sstauS.str();

   		Info << "Reading field tauS" << wtauS << endl;
    	tauS.set
   		(
   			id,
    		new volScalarField
  			(
  				IOobject
  				(
  					wtauS,
  					runTime.timeName(),
  					mesh,
  					IOobject::NO_READ,
  					IOobject::AUTO_WRITE
  				),
  				mesh,
  				dimensionedScalar(wtauS, dimTime, 1.0e10)
  			)
    	);
    }

/*
    // 熱泳動速度が粒子サイズに依存しないと仮定
    Info<< "Creating field Uth\n" << endl;
    const volVectorField Uth
    (
        IOobject
        (
            "Uth",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -0.55*mu*fvc::grad(T)/(rho*T)
    );
*/
// 粒子速度（泳動＋ガス）
    Info<< "Creating field Up\n" << endl;
    const volVectorField Up
    (
        IOobject
        (
            "Up",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );

    Info<< "Calculating face flux field UphiP\n" << endl;
    const surfaceScalarField UphiP
    (
        IOobject
        (
            "UphiP",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(Up) & mesh.Sf()
     );



//  scalarList momentSet(nMoment, 0.); // セルループ用モーメントセット
//    particleDatas::CONST_.setdeltaDb(T);



    //////////////////////////////////////////


    Info<< "Creating field nIteration" << endl;
    volScalarField nIteration
    (
        IOobject
        (
            "nIteration",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("nIteration", dimless, 0)
    );


#endif /* QMOMFOAM_CREATEQMOM_H_ */
