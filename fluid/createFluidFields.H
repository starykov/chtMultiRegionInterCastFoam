    // Initialise fluid field pointer lists
    PtrList<twoPhaseMixtureThermo> mixtureFluid(fluidRegions.size());
    PtrList<interfaceProperties> interfaceFluid(fluidRegions.size());
    PtrList<volScalarField> rhoFluid(fluidRegions.size());
    PtrList<volScalarField> TFluid(fluidRegions.size());
    PtrList<volScalarField> alpha1Fluid(fluidRegions.size());
    PtrList<volScalarField> alpha2Fluid(fluidRegions.size());
    PtrList<volScalarField> rho1Fluid(fluidRegions.size());
    PtrList<volScalarField> rho2Fluid(fluidRegions.size());
    PtrList<volScalarField> psi1Fluid(fluidRegions.size());
    PtrList<volScalarField> psi2Fluid(fluidRegions.size());
    PtrList<dimensionedScalar> pMinFluid(fluidRegions.size());
    PtrList<volScalarField> dgdtFluid(fluidRegions.size());
    PtrList<graph> gammaGraphFluid(fluidRegions.size());
    PtrList<graph> rho1GraphFluid(fluidRegions.size());
    PtrList<volScalarField> gammaFluid(fluidRegions.size());
    PtrList<volScalarField> KInvFluid(fluidRegions.size());
PtrList<dimensionedScalar> dFluid(fluidRegions.size());
//PtrList<volScalarField> KInvFluid(fluidRegions.size());
	PtrList<dimensionedScalar> TsolFluid(fluidRegions.size());
	PtrList<dimensionedScalar> TliqFluid(fluidRegions.size());
	PtrList<dimensionedScalar> LFluid(fluidRegions.size());
	PtrList<volScalarField> LSourceFluid(fluidRegions.size());
	PtrList<volScalarField> dgammadtFluid(fluidRegions.size());
	PtrList<volScalarField> gammaOldFluid(fluidRegions.size());

    PtrList<volVectorField> UFluid(fluidRegions.size());
    PtrList<volVectorField> gammaNFluid(fluidRegions.size());
    PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
    PtrList<surfaceScalarField> rhoPhiFluid(fluidRegions.size());
    PtrList<surfaceScalarField> phiCpFluid(fluidRegions.size());
    PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
    PtrList<compressible::turbulenceModel> turbulenceFluid(fluidRegions.size());
    PtrList<volScalarField> p_rghFluid(fluidRegions.size());
    PtrList<volScalarField> ghFluid(fluidRegions.size());
    PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
    PtrList<radiation::radiationModel> radiation(fluidRegions.size());
    PtrList<volScalarField> KFluid(fluidRegions.size());
    PtrList<volScalarField> dpdtFluid(fluidRegions.size());

    List<scalar> initialMassFluid(fluidRegions.size());
    List<bool> frozenFlowFluid(fluidRegions.size(), false);

    PtrList<fv::IOoptionList> fluidFvOptions(fluidRegions.size());

    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;


        Info<< "    Adding to p_rghFluid\n" << endl;
        p_rghFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "p_rgh",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to UFluid\n" << endl;
        UFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to phiFluid\n" << endl;
        phiFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(UFluid[i])
                    & fluidRegions[i].Sf()
            )
        );

        Info<< "    Adding to thermoFluid\n" << endl;

        mixtureFluid.set
        (
            i,
            //twoPhaseMixtureThermo::New(fluidRegions[i]).ptr()
			new twoPhaseMixtureThermo(fluidRegions[i])
        );


    	Info<< "    Adding to alpha1Fluid\n" << endl;
    	alpha1Fluid.set
    	(
    			i,
/*    			new volScalarField
    			(
    					IOobject
    					(
    							"alpha1",
    							runTime.timeName(),
    							fluidRegions[i],
    							IOobject::READ_IF_PRESENT,
    							IOobject::AUTO_WRITE
    					),*/
    					mixtureFluid[i].alpha1()
    			//)
    	);


	Info<< "    Adding to alpha2Fluid\n" << endl;
	alpha2Fluid.set
	(
			i,
		/*	new volScalarField
			(
					IOobject
					(
							"alpha2",
							runTime.timeName(),
							fluidRegions[i],
							IOobject::READ_IF_PRESENT,
							IOobject::AUTO_WRITE
					),*/
					mixtureFluid[i].alpha2()
			//)
	);
	//alpha2Fluid[i] = scalar(1.0) - alpha1Fluid[i];


	Info<< "    Adding to rho1Fluid\n" << endl;
	rho1Fluid.set
	(
			i,
			new volScalarField
			(
					IOobject
					(
							"rho1",
							runTime.timeName(),
							fluidRegions[i],
							IOobject::NO_READ,
							IOobject::AUTO_WRITE
					),
					//mixtureFluid[i].thermo1().rho()
					fluidRegions[i],
					dimensionedScalar
									                (
									                    "rho1",
														dimMass/dimVolume,
									                    1000
									                )
			)
	);

	Info<< "    Adding to rho2Fluid\n" << endl;
	rho2Fluid.set
	(
			i,
		/*	new volScalarField
			(
					IOobject
					(
							"rho2",
							runTime.timeName(),
							fluidRegions[i],
							IOobject::NO_READ,
							IOobject::AUTO_WRITE
					),*/
					mixtureFluid[i].thermo2().rho()
			//)
	);


	        Info<< "    Adding to TFluid\n" << endl;
	        TFluid.set
	        (
	            i,
	            /*
				new volScalarField
	            (
	                IOobject
	                (
	                    "T",
	                    runTime.timeName(),
	                    fluidRegions[i],
	                    IOobject::MUST_READ,
	                    IOobject::AUTO_WRITE
	                ),*/
	                mixtureFluid[i].T()
	            //)
	        );



    	        Info<< "    Adding to psi1Fluid\n" << endl;
    	        psi1Fluid.set
    	        (
    	            i,
    	            new volScalarField("psi1", mixtureFluid[i].thermo1().psi())
    	        );

    	        Info<< "    Adding to psi2Fluid\n" << endl;
    	            	        psi2Fluid.set
    	            	        (
    	            	            i,
    	            	            new volScalarField("psi2", mixtureFluid[i].thermo1().psi())
    	            	        );

        Info<< "    Adding to rhoFluid\n" << endl;
        rhoFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
				fluidRegions[i],
									dimensionedScalar
													                (
													                    "rho",
																		dimMass/dimVolume,
													                    1000
													                )
				//alpha1Fluid[i]*rho1Fluid[i] + alpha2Fluid[i]*rho2Fluid[i]
            )
        );

    	// Heat capacities
    	pMinFluid.set
    	(
    			i,
    			new dimensionedScalar(mixtureFluid[i].lookup("pMin"))
    	);
    	// Dendrite arm spacing
    	dFluid.set
    	(
    			i,
    			new dimensionedScalar(mixtureFluid[i].lookup("d"))
    	);

    	// Dendrite arm spacing
    	TsolFluid.set
    	(
    			i,
    			new dimensionedScalar(mixtureFluid[i].lookup("Tsol"))
    	);

    	// Dendrite arm spacing
    	TliqFluid.set
    	(
    			i,
    			new dimensionedScalar(mixtureFluid[i].lookup("Tliq"))
    	);

    	// Dendrite arm spacing
    	LFluid.set
    	(
    			i,
    			new dimensionedScalar(mixtureFluid[i].lookup("L"))
    	);


        Info<< "    Adding to gFluid\n" << endl;
        gFluid.set
        (
            i,
            new uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );


        Info<< "    Adding to ghFluid\n" << endl;
        ghFluid.set
        (
            i,
            new volScalarField("gh", gFluid[i] & fluidRegions[i].C())
        );

        Info<< "    Adding to ghfFluid\n" << endl;
        ghfFluid.set
        (
            i,
            new surfaceScalarField("ghf", gFluid[i] & fluidRegions[i].Cf())
        );

    	Info<< "    Adding to rhoPhiFluid\n" << endl;
    	rhoPhiFluid.set
    	(
    			i,
    			new surfaceScalarField
    			(
    					IOobject
    					(
    							"rhoPhi",
    							runTime.timeName(),
    							fluidRegions[i],
    							IOobject::READ_IF_PRESENT,
    							IOobject::AUTO_WRITE
    					),
    					//linearInterpolate(rhoFluid[i]*UFluid[i]) & fluidRegions[i].Sf()
						fvc::interpolate(rhoFluid[i])*phiFluid[i]
    			)
    	);



        Info<< "    Adding to dgdtFluid\n" << endl;
            	        dgdtFluid.set
            	        (
            	            i,
            	            new volScalarField(
            	            		IOobject
            	            		                (
            	            		                    "dgdt",
            	            		                    runTime.timeName(),
            	            		                    fluidRegions[i],
            	            		                    IOobject::NO_READ,
            	            		                    IOobject::AUTO_WRITE
            	            		                ),

            	            		pos(alpha2Fluid[i])*fvc::div(phiFluid[i])/max(alpha2Fluid[i], scalar(0.0001))
									)
            	        );




            	        Info<< "    Adding to interface\n" << endl;
            	        interfaceFluid.set
            	        (
            	            i,
							new interfaceProperties
            	            (
            	                alpha1Fluid[i],
            	                UFluid[i],
            	                mixtureFluid[i]
            	            )
            	        );


        Info<< "    Adding to turbulence\n" << endl;
        turbulenceFluid.set
        (
            i,
            compressible::turbulenceModel::New
            (
                rhoFluid[i],
                UFluid[i],
                rhoPhiFluid[i],
                mixtureFluid[i]
            ).ptr()
        );





        // Force p_rgh to be consistent with p
//        p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];

        radiation.set
        (
            i,
            radiation::radiationModel::New(mixtureFluid[i].T())
        );

        initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

        Info<< "    Adding to KFluid\n" << endl;
        KFluid.set
        (
            i,
            new volScalarField
            (
                "K",
                0.5*magSqr(UFluid[i])
            )
        );


        Info<< "    Adding to dpdtFluid\n" << endl;
        dpdtFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "dpdt",
                    runTime.timeName(),
                    fluidRegions[i]
                ),
                fluidRegions[i],
                dimensionedScalar
                (
                    "dpdt",
                    mixtureFluid[i].p().dimensions()/dimTime,
                    0
                )
            )
        );

	    IFstream file_gamma(runTime.path()/runTime.constant()/fluidRegions[i].name()/"gamma");
	    //Read file
	    gammaGraphFluid.set(
	    		i,
				new graph
	    (
	        "gamma_data_file",
	        "T_data",
	        "gamma_data",
	        file_gamma
	    ));

	    IFstream file_rho1(runTime.path()/runTime.constant()/fluidRegions[i].name()/"rho1");
	    //Read file
	    rho1GraphFluid.set(
	    		i,
				new graph
	    (
	        "rho1_data_file",
	        "T_data",
	        "rho1_data",
	        file_rho1
	    ));

        Info<< "    Adding to gammaFluid\n" << endl;
        gammaFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "gamma",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
				dimensionedScalar
				                (
				                    "gamma",
				                    dimless,
				                    0
				                )
            )
        );
        gammaFluid[i].oldTime();


        Info<< "    Adding to gammaFluid\n" << endl;
        LSourceFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "LSource",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
				dimensionedScalar
								                (
								                    "LSource",
													dimensionSet (1,-3,-1,1,0,0,0),
								                    0
								                )
            )
        );


        Info<< "    Adding to gammaFluid\n" << endl;
        dgammadtFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "dgammadt",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
				dimensionedScalar
								                (
								                    "dgammadt",
													dimensionSet (0,0,-1,0,0,0,0),
								                    0
								                )
            )
        );

        Info<< "    Adding to gammaFluid\n" << endl;
                gammaOldFluid.set
                (
                    i,
                    new volScalarField
                    (
                        IOobject
                        (
                            "gammaOld",
                            runTime.timeName(),
                            fluidRegions[i],
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        fluidRegions[i],
        				dimensionedScalar
        								                (
        								                    "gammaOld",
        													dimless,
        								                    0
        								                )
                    )
                );

        Info<< "    Adding to gammaFluid\n" << endl;
        KInvFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "KInv",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
				dimensionedScalar
				                (
				                    "KInv",
									dimless/sqr(dimLength),
				                    0
				                )
            )
        );



    	Info<< "    Adding to rho1Fluid\n" << endl;
    	gammaNFluid.set
    	(
    			i,
    			new volVectorField
    			(
    					IOobject
    					(
    							"gammaN",
    							runTime.timeName(),
    							fluidRegions[i],
    							IOobject::NO_READ,
    							IOobject::AUTO_WRITE
    					),
    					//mixtureFluid[i].thermo1().rho()
    					fluidRegions[i],
    					vector(0,0,0)
    			)
    	);
    	gammaNFluid[i]=fvc::grad(gammaFluid[i])*dimensionedScalar("one", dimLength, 1.0);
        gammaNFluid[i].dimensions().reset(dimless);
        gammaNFluid[i]=-1.0*pos(TFluid[i]-TsolFluid[i])*pos(TliqFluid[i]-TFluid[i])*gammaNFluid[i]/(mag(gammaNFluid[i])+VSMALL);


        const dictionary& pimpleDict =
            fluidRegions[i].solutionDict().subDict("PIMPLE");
        pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

        Info<< "    Adding fvOptions\n" << endl;
        fluidFvOptions.set
        (
            i,
            new fv::IOoptionList(fluidRegions[i])
        );
    }
