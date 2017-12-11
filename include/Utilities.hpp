//
// Created by matthias on 28.04.17.
//

#ifndef FIELDCALIBRATION_UTILITIES_H
#define FIELDCALIBRATION_UTILITIES_H

#endif //FIELDCALIBRATION_UTILITIES_H

#include "Laser.hpp"
#include "LaserTrack.hpp"


std::vector<Laser> ReachedExitPoint(const Laser&, float);

std::vector<Laser> SplitTrackSet(const Laser&, unsigned int);

Laser MergeLaser(const Laser &LaserA, const Laser &LaserB);

std::vector<Laser> InterlacedIterTrackSamples(const Laser &LaserSet);

LaserTrack Anode(TPCVolumeHandler &TPCVolume);