//
// Created by Yifan Chen on 22.01.18.
//

#ifndef FIELDCALIBRATION_VOXELMESH_HPP
#define FIELDCALIBRATION_VOXELMESH_HPP

#include "../include/ThreeVector.hpp"
#include "../include/LaserTrack.hpp"


std::vector<std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>>>
MeshVoxel(const std::vector<LaserTrack> &LaserTrackSet);





#endif //FIELDCALIBRATION_VOXELMESH_HPP
