dimensionedScalar dt("dt", dimTime, runTime.deltaT().value()/100);
dimensionedScalar localT("localT", dimTime, 0);
label nItt = 0;
bool timeComplete = false;
// U1 = U;
// U2 = U;
// thermo1.he() += magSqr(U)/2.0;
// thermo2.he() += magSqr(U)/2.0;

phase1.decode();
phase2.decode();
while (!timeComplete)
{
    tmp<volScalarField> Kd(fluid.Kd());

    volVectorField alphaRhoU1Old(alphaRhoU1);
    volVectorField alphaRhoU2Old(alphaRhoU2);
    volScalarField alphaRhoE1Old(alphaRhoE1);
    volScalarField alphaRhoE2Old(alphaRhoE2);

    //- 1st predictor
    phase1.encode();
    phase2.encode();

    volVectorField UK1(dt*Kd()*(U2 - U1));
    volScalarField EK1(dt*Kd()*(U & (U2 - U1)));

    alphaRhoU1 += 0.5*UK1;
    alphaRhoU2 -= 0.5*UK1;
    alphaRhoE1 += 0.5*EK1;
    alphaRhoE2 -= 0.5*EK1;

    phase1.decode();
    phase2.decode();
    p1 = neg(p1)*p2 + pos(p1)*p1;
    p2 = neg(p2)*p1 + pos(p2)*p2;

    U = fluid.U();

    //- 2nd predictor
    phase1.encode();
    phase2.encode();

    volVectorField UK2(dt*Kd()*(U2 - U1));
    volScalarField EK2(dt*Kd()*(U & (U2 - U1)));

    alphaRhoU1 += 0.5*UK2;
    alphaRhoU2 -= 0.5*UK2;
    alphaRhoE1 += 0.5*EK2;
    alphaRhoE2 -= 0.5*EK2;

    phase1.decode();
    phase2.decode();
    p1 = neg(p1)*p2 + pos(p1)*p1;
    p2 = neg(p2)*p1 + pos(p2)*p2;

    U = fluid.U();

    //- 3rd predictor
    phase1.encode();
    phase2.encode();

    volVectorField UK3(dt*Kd()*(U2 - U1));
    volScalarField EK3(dt*Kd()*(U & (U2 - U1)));

    alphaRhoU1 = alphaRhoU1Old + UK3;
    alphaRhoU2 = alphaRhoU2Old - UK3;
    alphaRhoE1 = alphaRhoE1Old + EK3;
    alphaRhoE2 = alphaRhoE2Old - EK3;

    phase1.decode();
    phase2.decode();
    p1 = neg(p1)*p2 + pos(p1)*p1;
    p2 = neg(p2)*p1 + pos(p2)*p2;

    volVectorField U1New(U1);
    volVectorField U2New(U2);

    U = fluid.U();

    //- Corrector
    phase1.encode();
    phase2.encode();

    volVectorField UK4
    (
        (UK1 + 2.0*UK2 + 2.0*UK3 + dt*Kd()*(U2 - U1))/6.0
    );
    volScalarField EK4
    (
        (EK1 + 2.0*EK2 + 2.0*EK3 + dt*Kd()*(U & (U2 - U1)))/6.0
    );

    alphaRhoU1 = alphaRhoU1Old + UK4;
    alphaRhoU2 = alphaRhoU2Old - UK4;
    alphaRhoE1 = alphaRhoE1Old + EK4;
    alphaRhoE2 = alphaRhoE2Old - EK4;

    phase1.decode();
    phase2.decode();
    p1 = neg(p1)*p2 + pos(p1)*p1;
    p2 = neg(p2)*p1 + pos(p2)*p2;

    U = fluid.U();

    scalar error = 0.0;
    error +=
        max
        (
            Foam::sqr(mag(U1New - U1))
           /(
                dimensionedScalar("ATol", dimVelocity, 1e-4)
              + 1e-6*mag(U1)
            )
        ).value();
    error +=
        max
        (
            Foam::sqr(mag(U2New - U2))
           /(
                dimensionedScalar("ATol", dimVelocity, 1e-8)
              + 1e-6*mag(U2)
            )
        ).value();
    error = Foam::sqrt(error);
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

Info<< "number of velocity relaxation iterations: " << nItt << endl;