//
// Created by Yifan Chen on 22.01.18.
//

#ifndef FIELDCALIBRATION_VOXELMESH_HPP
#define FIELDCALIBRATION_VOXELMESH_HPP

#include "../include/ThreeVector.hpp"
#include "../include/LaserTrack.hpp"


std::vector<std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>>>
MeshVoxel(const std::vector<LaserTrack> &LaserTrackSet,const TPCVolumeHandler &TPC);

std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>>
AveragebyDistance(std::vector<std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>>> &VoxelMesh,
                  const TPCVolumeHandler &TPC);





#endif //FIELDCALIBRATION_VOXELMESH_HPP
