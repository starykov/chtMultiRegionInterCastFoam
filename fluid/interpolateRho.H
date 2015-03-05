//gamma
  //Interpolation internal field
  gamma.field() = interpolateXY
  (
    T.field(),
    gammaGraph.x(),
	gammaGraph.y()
  )*alpha1;

  //Interpolation boundary field
  forAll(gamma.boundaryField(), patchi)
  {
	  gamma.boundaryField()[patchi] = interpolateXY
    (
      T.boundaryField()[patchi],
	  gammaGraph.x(),
	  gammaGraph.y()
    )*alpha1;
  }


  //rho1
    //Interpolation internal field
    rho1.field() = interpolateXY
    (
      T.field(),
      rho1Graph.x(),
	  rho1Graph.y()
    );

    //Interpolation boundary field
    forAll(rho1.boundaryField(), patchi)
    {
  	  rho1.boundaryField()[patchi] = interpolateXY
      (
        T.boundaryField()[patchi],
  	  rho1Graph.x(),
  	  rho1Graph.y()
      );
    }

