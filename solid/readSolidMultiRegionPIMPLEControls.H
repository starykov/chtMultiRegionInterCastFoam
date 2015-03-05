    const dictionary& pimples = mesh.solutionDict().subDict("PIMPLE");

    int nNonOrthCorr =
        pimples.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
