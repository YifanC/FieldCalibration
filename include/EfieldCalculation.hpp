//
// Created by Yifan Chen on 08.10.17.
//

#include "TROOT.h"
#include "TFile.h"
#include "TH3.h"
#include "../include/ThreeVector.hpp"
#include "../include/TPCVolumeHandler.hpp"
#include "../include/DriftVelocity.hpp"

#ifndef FIELDCALIBRATION_EFIELDCALCULATION_HPP
#define FIELDCALIBRATION_EFIELDCALCULATION_HPP

#endif //FIELDCALIBRATION_EFIELDCALCULATION_HPP

std::pair<std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>>
Efield(TPCVolumeHandler &TPCVolume, float cryoTemp, float E0, float v0, const char *root_name);

std::pair<std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>>
EfieldvecMap(TPCVolumeHandler &TPCVolume, float cryoTemp, float E0, float v0, std::vector<ThreeVector<float>> DMapTT);

std::tuple<std::vector<float >, std::vector<float >, std::vector<float >,
        std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>>
EfieldXYZwithBoundary(TPCVolumeHandler &TPCVolume, float cryoTemp, float E0, float v0, const char *root_name);

float EdgeEx(std::vector<ThreeVector<float>> &Efield, ThreeVector<unsigned long> Resolution, ThreeVector<float> Unit, float E0, ThreeVector<unsigned long> Coord);