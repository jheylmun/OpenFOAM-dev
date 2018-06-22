phaseModel& phase1 = fluid.phase1();
phaseModel& phase2 = fluid.phase2();

const volScalarField& alpha1 = phase1;
const volScalarField& alpha2 = phase2;

volVectorField& U1 = phase1.URef();
surfaceScalarField& phi1 = phase1.phiRef();
const surfaceScalarField& alphaPhi1 = phase1.alphaPhi();

volVectorField& U2 = phase2.URef();
surfaceScalarField& phi2 = phase2.phiRef();
const surfaceScalarField& alphaPhi2 = phase2.alphaPhi();

surfaceScalarField& phi = fluid.phi();

basicThermo& thermo1 = phase1.thermoRef();
basicThermo& thermo2 = phase2.thermoRef();

volScalarField& rho1 = phase1.rhoRef();
const volScalarField& psi1 = phase1.psi();

volScalarField& rho2 = phase2.rhoRef();
const volScalarField& psi2 = phase2.psi();

const IOMRFZoneList& MRF = fluid.MRF();
fv::options& fvOptions = fluid.fvOptions();
