//
// Created by Yifan Chen on 22.01.18.
//

#include "../include/WeightAverage.hpp"
#include "../include/Interpolation3D.hpp"
#include "../include/Laser.hpp"
#include "../include/ThreeVector.hpp"
#include <sstream>
#include <iostream>
#include <limits>
#include <vector>
#include <array>
#include <stdlib.h>


std::vector<std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>>>
MeshVoxel(const std::vector<LaserTrack> &LaserTrackSet,const TPCVolumeHandler &TPC){

    ThreeVector<unsigned long> Resolution = TPC.GetDetectorResolution();
    ThreeVector<float> MinimumCoord = TPC.GetMapMinimum();
    ThreeVector<float> MaximumCoord = TPC.GetMapMaximum();
    ThreeVector<float> Unit = {TPC.GetDetectorSize()[0] / (Resolution[0] - 1),
                               TPC.GetDetectorSize()[1] / (Resolution[1] - 1),
                               TPC.GetDetectorSize()[2] / (Resolution[2] - 1)};
    std::vector<std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>>> Mesh;
    Mesh.resize(Resolution[0]*Resolution[1]*Resolution[2]);

    // Loop over data points (tracks) of the whole sample
    for (unsigned long track = 0; track < LaserTrackSet.size(); track++) {
        // Loop over data points (samples) of each track
        for (unsigned long sample = 0; sample < LaserTrackSet[track].GetNumberOfSamples(); sample++) {

            ThreeVector<float> Position = LaserTrackSet[track].GetSamplePosition(sample);
            ThreeVector<float> Disp = LaserTrackSet[track].GetDisplacement(sample);

            // x,y,z,dx,dy,dz
//            std::pair<ThreeVector<float >, ThreeVector<float>>> SampleInfo = std::make_pair(Position,Disp);
            auto SampleInfo = std::make_pair(Position,Disp);

            int xbin = (Position[0] - MinimumCoord[0] - Unit[0]*0.5) / Unit[0];
            int ybin = (Position[1] - MinimumCoord[1] - Unit[1]*0.5) / Unit[1];
            int zbin = (Position[2] - MinimumCoord[2] - Unit[2]*0.5) / Unit[2];

            // The order of the vector given by "zbin + (Resolution[2] * (ybin + Resolution[1] * xbin))"
            // should be consistent for writing and reading
            Mesh[zbin + (Resolution[2] * (ybin + Resolution[1] * xbin))].push_back(SampleInfo);

        } // end sample loop
    } // end track loop

    return Mesh;

};


std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>>
AveragebyDistance(std::vector<std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>>> &VoxelMesh,
                  const TPCVolumeHandler &TPC){

    ThreeVector<unsigned long> Resolution = TPC.GetDetectorResolution();
    ThreeVector<float> Unit = {TPC.GetDetectorSize()[0] / (Resolution[0] - 1),
                               TPC.GetDetectorSize()[1] / (Resolution[1] - 1),
                               TPC.GetDetectorSize()[2] / (Resolution[2] - 1)};
    ThreeVector<float> MinimumCoord = TPC.GetMapMinimum();

    //Initialize the map of displacement and its error(deviation)
    std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>> AverageGrid;

    for(int i = 0; i<VoxelMesh.size(); i++){

        // zbin + (Resolution[2] * (ybin + Resolution[1] * xbin))
        int xbin = div(i,Resolution[1]*Resolution[2]).quot;
        int ybin = div(div(i,Resolution[1]*Resolution[2]).rem, Resolution[2]).quot;
        int zbin = div(div(i,Resolution[1]*Resolution[2]).rem, Resolution[2]).rem;

        ThreeVector<float> Grid = {MinimumCoord[0] + Unit[0] * xbin,
                                   MinimumCoord[1] + Unit[1] * ybin,
                                   MinimumCoord[2] + Unit[2] * zbin};

        // start weighed mean calculation
        std::vector<float> vecw;
        std::vector<ThreeVector<float>> vecDisp;

        float sumw = 0;
        ThreeVector<float> sumDispw = {0, 0, 0};

        for(int j = 0; j <VoxelMesh[i].size(); j++){

            ThreeVector<float> SamplePosition = VoxelMesh[i][j].first;
            ThreeVector<float> SampleDistortion = VoxelMesh[i][j].second;

            // Take 1/R as weight
            // Shall we try 1/R2 as well?
            float w = 1 / sqrtf((Grid[0]-SamplePosition[0])*(Grid[0]-SamplePosition[0])
                            +(Grid[1]-SamplePosition[1])*(Grid[1]-SamplePosition[1])
                            +(Grid[2]-SamplePosition[2])*(Grid[2]-SamplePosition[2]));


            vecDisp.push_back(SampleDistortion);

            vecw.push_back(w);

            //Normalization factor
            sumw += w;

            sumDispw += SampleDistortion * w;
        }

        ThreeVector<float> Average = sumDispw/sumw;
        //End of weighed mean calculation

        //Start for weighed deviation
        ThreeVector<float> sigma = {0, 0, 0};

        for(int j = 0; j < VoxelMesh[i].size(); j++){
            for(int k = 0; k < 3; k++){
                sigma[k] += vecw[j]/sumw * (vecDisp[j][k]-Average[k])* (vecDisp[j][k]-Average[k]);
            }
        }

        ThreeVector<float> Deviation = {sqrtf(sigma[0]), sqrtf(sigma[1]), sqrtf(sigma[2])};
        //End for weighed deviation

        std::pair<ThreeVector<float >, ThreeVector<float>> DistortionError = std::make_pair(Average, Deviation);

        AverageGrid.push_back(DistortionError);

    }

    return AverageGrid;

};
