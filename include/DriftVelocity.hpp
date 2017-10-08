//
// Created by Yifan Chen on 28.09.17.
//

#ifndef FIELDCALIBRATION_DRIFTVELOCITY_HPP
#define FIELDCALIBRATION_DRIFTVELOCITY_HPP

#endif //FIELDCALIBRATION_DRIFTVELOCITY_HPP

float ElectronDriftVelocity(float T, float E);

float elecDriftVelHelper(float T, float E, int whichFit);

float searchE(float v_drift, float cryoTemp, float E0);