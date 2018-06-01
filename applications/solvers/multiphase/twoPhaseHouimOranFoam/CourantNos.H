//     #include "CourantNo.H"

{
    CoNum = 0.5*gMax
    (
        fvc::surfaceSum
        (
            max
            (
                mag(phi1)
              + fvc::interpolate
                (
                    phase1.c()
                )*mesh.magSf(),
                 mag(phi2)
              + fvc::interpolate
                (
                    phase2.c()
                )*mesh.magSf()
            )
        )().primitiveField()/mesh.V().field()
    )*runTime.deltaTValue();

    Info<< "Max Ur Courant Number = " << CoNum << endl;
}
