//
// Created by Yifan Chen on 09.10.17.
//

#ifndef FIELDCALIBRATION_LASERDISPLACEMENTITERATION_HPP
#define FIELDCALIBRATION_LASERDISPLACEMENTITERATION_HPP

#endif //FIELDCALIBRATION_LASERDISPLACEMENTITERATION_HPP

#include "../include/Laser.hpp"

std::pair<Laser, Laser> DispLaserIteration(unsigned Nstep, Laser LaserSet1, Laser LaserSet2, bool CorrMapFlag);