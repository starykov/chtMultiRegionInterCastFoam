//gamma
  //Interpolation internal field
  gamma.field() = interpolateXY
  (
    T.field(),
    gammaGraph.x(),
	gammaGraph.y()
  );

  //Interpolation boundary field
  forAll(gamma.boundaryField(), patchi)
  {
	  gamma.boundaryField()[patchi] = interpolateXY
    (
      T.boundaryField()[patchi],
	  gammaGraph.x(),
	  gammaGraph.y()
    );
  }


