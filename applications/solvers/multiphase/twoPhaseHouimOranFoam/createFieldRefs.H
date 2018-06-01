volScalarField& rho1 = phase1.rho();
volScalarField& rho2 = phase2.rho();

volScalarField& p1 = phase1.p();
volScalarField& E1 = phase1.E();

volScalarField& p2 = phase2.p();
volScalarField& E2 = phase2.E();

rhoThermo& thermo1 = phase1.thermo();
volScalarField& he1 = thermo1.he();

rhoThermo& thermo2 = phase2.thermo();
volScalarField& he2 = thermo2.he();
