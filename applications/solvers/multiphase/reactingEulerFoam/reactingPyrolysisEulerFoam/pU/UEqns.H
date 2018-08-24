Info<< "Constructing momentum equations" << endl;

PtrList<fvVectorMatrix> UEqns(phases.size());


const phaseModel& solidPhase = phases[solidReactionPhases[0]];
const phaseModel& gasPhase = phases[gasReactionPhases[0]];

volScalarField sumRRs(solidPhase.dmdt());
volScalarField xi(pos0(sumRRs));

{
    autoPtr<phaseSystem::momentumTransferTable>
        momentumTransferPtr(fluid.momentumTransfer());

    phaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    forAll(fluid.movingPhases(), movingPhasei)
    {
        phaseModel& phase = fluid.movingPhases()[movingPhasei];

        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();
        volVectorField& U = phase.URef();

        UEqns.set
        (
            phase.index(),
            new fvVectorMatrix
            (
                phase.UEqn()
             ==
               *momentumTransfer[phase.name()]
              + fvOptions(alpha, rho, U)
            )
        );

        if (solidReactionPhases.found(phase.name()))
        {
            UEqns[phase.index()] +=
                fvm::Sp(sumRRs*xi, U)
              - sumRRs*(1.0 - xi)*gasPhase.U();
        }
        else if (gasReactionPhases.found(phase.name()))
        {
            UEqns[phase.index()] +=
                fvm::Sp(sumRRs*(1.0 - xi), U)
              - sumRRs*xi*solidPhase.U();
        }


        UEqns[phase.index()].relax();
        fvOptions.constrain(UEqns[phase.index()]);
        fvOptions.correct(U);
    }
}
