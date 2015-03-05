    //const scalar little = 1e-3;

	fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + turbulence.divDevRhoReff(U)
    );

	UEqn += fvm::Sp(rho*mixture.thermo1().nu()*KInv, U);
	

    UEqn.relax();

    if (momentumPredictor)
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    interface.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        K = 0.5*magSqr(U);
    }
