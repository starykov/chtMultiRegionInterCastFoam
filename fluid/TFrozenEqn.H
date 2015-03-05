{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
//      + fvm::div(rhoPhi, T)
      - fvm::laplacian(mixture.alphaEff(turbulence.mut()), T)
/*    + (
            fvc::div(fvc::absolute(phi, U), p)
          + fvc::ddt(rho, K) + fvc::div(rhoPhi, K)
        )
       *(
           alpha1/mixture.thermo1().Cv()
         + alpha2/mixture.thermo2().Cv()
        )
*/
    );

    TEqn.relax();
    TEqn.solve();

    //mixture.T()=T;

    mixture.correct();

    //mixture.thermo1().he = mixture.thermo1().he
    //mixture.thermo1().correct();
    //mixture.thermo2().correct();

    Info<< "min(T) " << min(T).value() << " max(T) " << max(T).value() << endl;
    Info<< "min(T(alloy)) " << min(T*alpha1).value() << " max(T(alloy)) " << max(T*alpha1).value() << endl;
}
