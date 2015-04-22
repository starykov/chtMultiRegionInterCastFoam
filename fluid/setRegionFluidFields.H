

	fvMesh& mesh = fluidRegions[i];
	Info << "region number " << i <<endl;


	twoPhaseMixtureThermo& mixture = mixtureFluid[i];
    //thermo.validate(args.executable(), "h" ,"e");


    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    //volScalarField& T = TFluid[i];


    volScalarField& alpha1 = alpha1Fluid[i];
    volScalarField& alpha2 = alpha2Fluid[i];


    compressible::turbulenceModel& turbulence = turbulenceFluid[i];
    volScalarField& K = KFluid[i];
    volScalarField& dpdt = dpdtFluid[i];

    volScalarField& dgdt = dgdtFluid[i];

    volScalarField& p = mixture.p();
    const volScalarField& psi = mixture.psi();

    volScalarField& T = mixture.T();
    //volScalarField& rho1 = rho1Fluid[i];//
    volScalarField& rho1 = mixture.thermo1().rho();//rho1Fluid[i];
    const volScalarField& psi1 = psi1Fluid[i];
    volScalarField& rho2 = mixture.thermo2().rho();//rho2Fluid[i];
    const volScalarField& psi2 = psi2Fluid[i];

    graph& rho1Graph = rho1GraphFluid[i];
    //rho1.write();
    graph& gammaGraph = gammaGraphFluid[i];

    volScalarField& gamma = gammaFluid[i];

 	#include "interpolateGamma.H"
	//#include "interpolateRho.H"

    //rho1 = mixture.thermo1().rho();

    volScalarField& rho = rhoFluid[i];
    rho = alpha1*rho1 + alpha2*rho2;
    surfaceScalarField& rhoPhi = rhoPhiFluid[i];


    volScalarField& KInv = KInvFluid[i];
	dimensionedScalar& d = dFluid[i];
	dimensionedScalar& Tsol = TsolFluid[i];
	dimensionedScalar& Tliq = TliqFluid[i];
	dimensionedScalar& L = LFluid[i];

	volVectorField& gammaN = gammaNFluid[i];
	volScalarField& LSource = LSourceFluid[i];
	volScalarField& dgammadt = dgammadtFluid[i];
	volScalarField& gammaOld = gammaOldFluid[i];
	volScalarField& dgamma = dgammaFluid[i];
	volScalarField& v = vFluid[i];


    //rhoPhi.write();

    dimensionedScalar& pMin = pMinFluid[i];
    interfaceProperties& interface = interfaceFluid[i];
    //interfaceProperties interface(alpha1, U, mixture);

    volScalarField& p_rgh = p_rghFluid[i];
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    radiation::radiationModel& rad = radiation[i];

    fv::IOoptionList& fvOptions = fluidFvOptions[i];



    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    bool& frozenFlow = frozenFlowFluid[i];
