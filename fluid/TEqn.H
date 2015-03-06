    int iter=0;
    word solver, field;
    double ini=0.0, fin=0.0;
    int niter=0.0;
    double diff1=0.0;
    double diff2=0.0;
    double diff3=0.0;
    volScalarField gammaO("gammaO",gamma);

    solverPerformance::debug=0; // do not write solver specification out again

    gamma.storePrevIter();

    do
      {

    	gammaOld = gamma.prevIter();
    	 gamma.storePrevIter();


    	    iter++;
    	//if (gamma == gammaOld) {Info << "gamma is equal to gammaold" << endl;}
    	    dgammadt = (gamma-gammaOld)/runTime.deltaT();
    LSource = L*rho/mixture.Cp()*(dgammadt);



    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - fvm::laplacian(mixture.alphaEff(turbulence.mut()), T)

    + (
            fvc::div(fvc::absolute(phi, U), p)
          + fvc::ddt(rho, K) + fvc::div(rhoPhi, K)
        )
       *(
           alpha1/mixture.thermo1().Cv()
         + alpha2/mixture.thermo2().Cv()
        )
		== LSource//L*rho/mixture.Cp()*(fvc::ddt(gamma))
    );

    TEqn.relax();
    //lduMatrix::
    Foam::solverPerformance sp=TEqn.solve();
      ini=sp.initialResidual();
      fin=sp.finalResidual();
      niter+=sp.nIterations();

    #include "interpolateGamma.H"
    gamma.relax();


    diff1 = Foam::gMax(mag(gamma.internalField()-gamma.prevIter().internalField()));




    if (iter==1) {diff3=ini; solver=sp.solverName(); field=sp.fieldName(); diff2=diff1;}

  }while ((iter < iterMax) && (diff1 > convergence));


    //TEqn.relax();
    //TEqn.solve();

    //mixture.T()=T;

    mixture.correct();

    solverPerformance::debug=1; // write out again solver news

      Info<< solver << ": " << " Solving for " << field << ", Initial residual =  " << diff3 << ", Final residual = " << fin << ", No Iterations "<< niter << " SUM" << endl;

      Info<< "SUBCYCLE:  Solving for fs, Initial residual =  " << diff2 << ", Final residual = " << diff1 << ", No Iterations "<< iter << endl;

      if (iter==iterMax) {Info<< "### WARNING: fs subcycle reached maximum iteration of "<< iterMax << " ###" << nl << endl;}

      //rhok = 1.0 - beta*(T - Tliq);
      //v = pos(T-Tsol)*pos(Tliq-T)*mag(fvc::ddt(fs))/(mag(fvc::grad(fs))+dimensionedScalar("SMALL", dimless/dimLength, SMALL));

      // normal vector for the solid fraction contours
      gammaN=fvc::grad(gamma)*dimensionedScalar("one", dimLength, 1.0);
      gammaN.dimensions().reset(dimless);
      gammaN=-1.0*pos(T-Tsol)*pos(Tliq-T)*gammaN/(mag(gammaN)+SMALL);

    Info<< "min(T) " << min(T).value() << " max(T) " << max(T).value() << endl;
    Info<< "min(T(alloy)) " << min(T*alpha1).value() << " max(T(alloy)) " << max(T*alpha1).value() << endl;

