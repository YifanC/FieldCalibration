//
// Created by Yifan Chen on 28.09.17.
//
#include <sstream>
#include <iostream>
#include <limits>
#include "../include/DriftVelocity.hpp"
#include <limits>
#include <cmath>

//T in K, E in kV/cm, Drift velocity in mm/us
float ElectronDriftVelocity(float T, float E) {
    float results;

    float xFit = 0.0938163 - 0.0052563 * (T - 87.302) - 0.000146981 * std::pow(T - 87.302, 2);
    float muFit = 5.183987 + 0.01447761 * (T - 87.302) - 0.0034972 * std::pow(T - 87.302, 2) -
                  0.0005162374 * std::pow(T - 87.302, 3);

    if (E < xFit)
        results = E * muFit;
    else if (E <= 0.619)
        results = elecDriftVelHelper(T, E, 2);
    else if (E >= 0.699)
        results = elecDriftVelHelper(T, E, 1);
    else
        results = (E - 0.619) / 0.08 * elecDriftVelHelper(T, E, 1) + (0.699 - E) / 0.08 * elecDriftVelHelper(T, E, 2);

    return results;
}

//T in K, E in kV/cm
float elecDriftVelHelper(float T, float E, int whichFit) {
    float p1 = 0.0;
    float p2 = 0.0;
    float p3 = 0.0;
    float p4 = 0.0;
    float p5 = 0.0;
    float p6 = 0.0;
    float t0 = 0.0;

    if (whichFit == 1)  //WalkowiakParameterSet
    {
        p1 = -0.01481;
        p2 = -0.0075;
        p3 = 0.141;
        p4 = 12.4;
        p5 = 1.627;
        p6 = 0.317;
        t0 = 90.371;
    } else if (whichFit == 2)  //IcarusFit
    {
        p1 = -0.04640231;
        p2 = 0.0171171;
        p3 = 1.881246;
        p4 = 0.9940772;
        p5 = 0.0117183;
        p6 = 4.202141;
        t0 = 105.7491;
    }

    float results;
    results = (1 + p1 * (T - t0)) * (p3 * E * log(1 + fabs(p4) / E) + p5 * std::pow(E, p6)) + p2 * (T - t0);

    return results;
}


// v_drift(E) is function of E (monotonically increasing)
// v_drift (mm/us) ; E (kV/cm)
// This function serves as "inverse function" for v_drift(E)
// It returns the magnitude of E field with given v_drift with precision XXXX
float searchE(float v_drift, float cryoTemp, float E0) {
    float Eresult = std::numeric_limits<float>::max();
    float Emin = E0 * (1 - 1); // The distortion of E field reaches 100%
    float Emax = E0 * (1 + 1); // The distortion of E field reaches 100%, the two side is not exactly the same

    // the values of HalfTolerance and Nmax are tricky, the precision of v may not go over 1E-3
    float HalfTolerance = 1E-3; // This Tolerence could/should have the same order as the error size from the fitting model
    const int Nmax = 20;
    int n = 0; // count on the dividing steps, in case of dead loop 2^20 ~ 1E6

    if (v_drift > ElectronDriftVelocity(cryoTemp, Emax) || v_drift < ElectronDriftVelocity(cryoTemp, Emin)) {
        //Something is very very very wrong
        std::cout
                << "The given drift velocity refers to a E field which is biased more than 100% E0 in the given condition."
                << std::endl;
    } else {
        // Bisection method
        while (Eresult != (Emin + Emax) * 0.5 && n < Nmax) {
            if (std::abs(v_drift - ElectronDriftVelocity(cryoTemp, (Emin + Emax) * 0.5)) < HalfTolerance) {
                Eresult = (Emin + Emax) * 0.5;
                break;
            } else if (v_drift > ElectronDriftVelocity(cryoTemp, (Emin + Emax) * 0.5) + HalfTolerance) {
                Emin = (Emin + Emax) * 0.5;
            }
                // This two sides is not exactly symmetric
            else if (v_drift < ElectronDriftVelocity(cryoTemp, (Emin + Emax) * 0.5) - HalfTolerance) {
                Emax = (Emin + Emax) * 0.5;
            }
            n++;
        }
    }

    if (Eresult != (Emin + Emax) * 0.5) {
        // Something is very very very wrong
        std::cout << "With drift velocity: " << v_drift
                  << "(mm/us), this function failed to find corresponding E field in the given condition. E: "
                  << 0.5 * (Emin + Emax)
                  << ", the v calculated is " << ElectronDriftVelocity(cryoTemp, 0.5 * (Emin + Emax))
                  << " and the nstep is " << n << std::endl;
    }

    return Eresult;
}
