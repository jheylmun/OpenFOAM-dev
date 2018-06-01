dimensionedScalar localT("localT", dimTime, 0);
label nItt = 0;
dimensionedScalar tauP("tauP", dimensionSet(0, 0, 1, 0, 0 , 0, 0), 5e-4);
dimensionedScalar dt
(
    "dt",
    dimTime,
    runTime.deltaT().value()/100
);
bool timeComplete = false;

while (!timeComplete)
{

    // Store old fields
    volScalarField alpha1Old(alpha1);
    volVectorField alphaRhoU1Old(alphaRhoU1);
    volVectorField alphaRhoU2Old(alphaRhoU2);
    volScalarField alphaRhoE1Old(alphaRhoE1);
    volScalarField alphaRhoE2Old(alphaRhoE2);

    // 1st Predictor
    phase1.encode();
    phase2.encode();

    volScalarField theta(alpha1*alpha2/(tauP*(p1 + p2)));
    volScalarField alphaK1(dt*theta*(p2 - p1));
    volScalarField alphaRhoEK1(dt*p*theta*(p2 - p1));

    alpha1 -= 0.5*alphaK1;
    alpha2 = 1.0 - alpha1;

    alphaRhoE1 += 0.5*alphaRhoEK1;
    alphaRhoE2 -= 0.5*alphaRhoEK1;

    phase1.decode();
    phase2.decode();
    p1 = neg(p1)*p2 + pos(p1)*p1;
    p2 = neg(p2)*p1 + pos(p2)*p2;
    p = fluid.p();

    // 2nd Predictor
    phase1.encode();
    phase2.encode();

    theta = alpha1*alpha2/(tauP*(p1 + p2));
    volScalarField alphaK2(dt*theta*(p2 - p1));
    volScalarField alphaRhoEK2(dt*p*theta*(p2 - p1));

    alpha1 -= 0.5*alphaK2;
    alpha2 = 1.0 - alpha1;

    alphaRhoE1 += 0.5*alphaRhoEK2;
    alphaRhoE2 -= 0.5*alphaRhoEK2;

    phase1.decode();
    phase2.decode();
    p1 = neg(p1)*p2 + pos(p1)*p1;
    p2 = neg(p2)*p1 + pos(p2)*p2;
    p = fluid.p();

    // 3rd Predictor
    phase1.encode();
    phase2.encode();

    theta = alpha1*alpha2/(tauP*(p1 + p2));
    volScalarField alphaK3(dt*theta*(p2 - p1));
    volScalarField alphaRhoEK3(dt*p*theta*(p2 - p1));

    alpha1 = alpha1Old - alphaK3;
    alpha2 = 1.0 - alpha1;

    alphaRhoE1 = alphaRhoE1Old + alphaRhoEK3;
    alphaRhoE2 = alphaRhoE2Old - alphaRhoEK3;

    phase1.decode();
    phase2.decode();
    p1 = neg(p1)*p2 + pos(p1)*p1;
    p2 = neg(p2)*p1 + pos(p2)*p2;
    p = fluid.p();

    volScalarField alphaNew(alpha1);
    volScalarField pNew(p);

    // Corrector
    phase1.encode();
    phase2.encode();

    theta = alpha1*alpha2/(tauP*(p1 + p2));
    volScalarField alphaK4
    (
        (
            alphaK1 + 2.0*alphaK2 + 2.0*alphaK3
          + dt*theta*(p2 - p1)
        )/6.0
    );
    volScalarField alphaRhoEK4
    (
        (
            alphaRhoEK1 + 2.0*alphaRhoEK2 + 2.0*alphaRhoEK3
          + dt*p*theta*(p2 - p1)
        )/6.0
    );

    alpha1 = alpha1Old - alphaK4;
    alpha2 = 1.0 - alpha1;

    alphaRhoE1 = alphaRhoE1Old + alphaRhoEK4;
    alphaRhoE2 = alphaRhoE2Old - alphaRhoEK4;

    phase1.decode();
    phase2.decode();
    p1 = neg(p1)*p2 + pos(p1)*p1;
    p2 = neg(p2)*p1 + pos(p2)*p2;
    p = fluid.p();

    scalar error =
        max
        (
            mag(pNew - p)
           /(
                dimensionedScalar("ATol", dimPressure, 1-4)
              + 1e-6*max(p, pNew)
            )
        ).value();
    error = max(error, SMALL);

    if (error < 1)
    {
        dt *= min(2, max(0.5, 0.9/Foam::pow(error, 1.0/3.0)));
        dimensionedScalar maxLocalDt =
            max
            (
                runTime.deltaT() - localT,
                dimensionedScalar("0", dimTime, 0.0)
            );
        dt = min(maxLocalDt, dt);
        if (dt.value() == 0.0)
        {
            timeComplete = true;
        }
        localT += dt;
    }
    else
    {
        dt *= min(1, max(0.5, 0.9/Foam::pow(error, 1.0/3.0)));

        alphaRhoU1 = alphaRhoU1Old;
        alphaRhoU2 = alphaRhoU2Old;
        alphaRhoE1 = alphaRhoE1Old;
        alphaRhoE2 = alphaRhoE2Old;
        phase1.decode();
        phase2.decode();
    }
    nItt++;
}
Info<< "number of pressure relaxation iterations: " << nItt << endl;