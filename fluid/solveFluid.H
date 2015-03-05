if (finalIter)
{
	mesh.data::add("finalIteration", true);
}

if (max(T*alpha1) < Tsol)
	{frozenFlow = true;
	Info << "flow is frozen, only temperature is calculated" << endl;
	}


const scalar little = 1e-3;
//mixture.thermo1().rho() = rho1;

KInv=sqr(gamma)/(pow3(scalar(1.0)-gamma)+little)/sqr(d);

if (frozenFlow)
{
	#include "TFrozenEqn.H"
}

else {

	#include "alphaEqnsSubCycle.H"

	interface.correct();

	mixture.alpha1()=alpha1;
	mixture.alpha2()=alpha2;

	solve(fvm::ddt(rho) + fvc::div(rhoPhi));

	#include "UEqn.H"
//    #include "EEqn.H"

	#include "TEqn.H"

	// --- PISO loop
	for (int corr=0; corr<nCorr; corr++)
	{
	#include "pEqn.H"
	}

	turbulence.correct();
}

//    rho = mixture.rho();

if (finalIter)
{
	mesh.data::remove("finalIteration");
}
