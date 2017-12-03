//
// Created by matthias on 28.04.17.
//

#ifndef FIELDCALIBRATION_UTILITIES_H
#define FIELDCALIBRATION_UTILITIES_H

#endif //FIELDCALIBRATION_UTILITIES_H

#include "Laser.hpp"

std::vector<Laser> ReachedExitPoint(const Laser&, float);

std::vector<Laser> SplitTrackSet(const Laser&, unsigned int);

std::vector<Laser> InterlacedIterTrackSamples(const Laser &LaserSet);