//
// Created by Yifan Chen on 22.01.18.
//

#include "../include/VoxelMesh.hpp"
#include "../include/Interpolation3D.hpp"
#include "../include/Laser.hpp"
#include "../include/ThreeVector.hpp"
#include <sstream>
#include <iostream>
#include <limits>
#include <vector>
#include <array>


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

/*
ThreeVector<float>
GridInterpolation(const std::vector<std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>>> Mesh,
                  ThreeVector<float> Location, const TPCVolumeHandler &TPC){
    ThreeVector<unsigned long> Resolution = TPC.GetDetectorResolution();
    ThreeVector<float> MinimumCoord = TPC.GetMapMinimum();
    ThreeVector<float> MaximumCoord = TPC.GetMapMaximum();
    ThreeVector<float> Unit = {TPC.GetDetectorSize()[0] / (Resolution[0] - 1),
                               TPC.GetDetectorSize()[1] / (Resolution[1] - 1),
                               TPC.GetDetectorSize()[2] / (Resolution[2] - 1)};
    int xbin = (Location[0] - MinimumCoord[0] - Unit[0]*0.5) / Unit[0];
    int ybin = (Location[1] - MinimumCoord[1] - Unit[1]*0.5) / Unit[1];
    int zbin = (Location[2] - MinimumCoord[2] - Unit[2]*0.5) / Unit[2];

    Location[0] = TPC.GetDetectorOffset()[0] + Unit[0] * xbin;
    Location[1] = TPC.GetDetectorOffset()[1] + Unit[1] * ybin;
    Location[2] = TPC.GetDetectorOffset()[2] + Unit[2] * zbin;


    float float_max = std::numeric_limits<float>::max();

    std::array<int, 4> VertexVoxel = {zbin + (Resolution[2] * ((ybin+1) + Resolution[1] * xbin)),
                                      (zbin+1) + (Resolution[2] * ((ybin-1) + Resolution[1] * xbin)),
                                      (zbin-1) + (Resolution[2] * ((ybin-1) + Resolution[1] * (xbin+1))),
                                      (zbin-1) + (Resolution[2] * ((ybin-1) + Resolution[1] * (xbin-1)))};

    //if none of the 4 bins is empty
    for(int m = 0; m < Mesh[VertexVoxel[0]].size(); m++){
        for(int n = 0; n < Mesh[VertexVoxel[1]].size(); n++){
            for(int l = 0; l < Mesh[VertexVoxel[2]].size(); l++){
                for(int k = 0; k < Mesh[VertexVoxel[3]].size(); k++){

                    // Initialize a displacement vector with zero
                    ThreeVector<float> InterpolatedDispl = {0.0, 0.0, 0.0};

                    // Initialize Barycentric coordinate system (it will have 4 dimensions)
                    std::vector<float> BaryCoord;

                    // Initialize matrix for Location transformation into barycentric coordinate system
                    Matrix3x3 TransMatrix = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
                    for (unsigned row = 0; row < 3; row++) {
                        // Loop over matrix columns
                        for (unsigned column = 0; column < 3; column++) {
                            // Fill transformation matrix elements
                            // x1,2,3 - x4, y1,2,3 - y4, z1,2,3 - z4
                            // 1,2,3,4 order is random
                            TransMatrix[row][column] =
                                    Mesh[VertexVoxel[row]][m][row] -
                                    Mesh[PointIndex.back().first].GetSamplePosition(PointIndex.back().second)[row];
                            TransMatrix[row][column] =
                                    LaserMeshSet[PointIndex[column].first].GetSamplePosition(PointIndex[column].second)[row] -
                                    LaserMeshSet[PointIndex.back().first].GetSamplePosition(PointIndex.back().second)[row];
                        }
                    }
                }
            }
        }
    }

    // Create a array which contains the info of all 4 vertices of a cell
    !!!!!!!!!std::array<std::pair<unsigned long, unsigned long>, 4> PointIndex;

    // Initialize a displacement vector with zero
    ThreeVector<float> InterpolatedDispl = {0.0, 0.0, 0.0};

    // Initialize Barycentric coordinate system (it will have 4 dimensions)
    std::vector<float> BaryCoord;

    // Find cell in the mesh where the point is located
    Delaunay::Cell_handle Cell = Mesh.locate(VectorToPoint(Location));

    // Loop over all four vertex points of the cell of interest
    for (unsigned vertex_no = 0; vertex_no < PointIndex.size(); vertex_no++) {
        // Get vertex info of the cell (track number, sample number)
        PointIndex[vertex_no] = Cell->vertex(vertex_no)->info();
    }

    // Initialize matrix for Location transformation into barycentric coordinate system
    Matrix3x3 TransMatrix = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    // Loop over matrix rows
    for (unsigned row = 0; row < 3; row++) {
        // Loop over matrix columns
        for (unsigned column = 0; column < 3; column++) {
            // Fill transformation matrix elements
            // x1,2,3 - x4, y1,2,3 - y4, z1,2,3 - z4
            TransMatrix[row][column] =
                    LaserMeshSet[PointIndex[column].first].GetSamplePosition(PointIndex[column].second)[row] -
                    LaserMeshSet[PointIndex.back().first].GetSamplePosition(PointIndex.back().second)[row];
        }
    }

    // Reuse Location and store its position relative to the last vertex of the cell it is contained in
    // after is step Location is (r-r4)
    Location -= LaserMeshSet[PointIndex.back().first].GetSamplePosition(PointIndex.back().second);

    // If the transformation matrix can be successfully inverted
    if (TransMatrix.Invert()) {
        // Use inverted matrix to fill the first three coordinates
        ThreeVector<float> BC = TransMatrix * Location;
        BaryCoord = BC.GetStdVector();

        // The sum of all barycentric coordinates has to be 1 by definition, use this to calculate the 4th coordinate
        BaryCoord.push_back(1 - BaryCoord[0] - BaryCoord[1] - BaryCoord[2]);
    }
    else // if the matrix can't be inverted
    {
        // Set displacement zero and end function immediately!
//        std::cout<<"The transition matrix for this D grid point is not invertable. "<<std::endl;
        InterpolatedDispl = {float_max, float_max, float_max};
//        if (Map) { InterpolatedDispl = {float_max, float_max, float_max}; }
//        else { InterpolatedDispl = {0, 0, 0}; }
//        InterpolatedDispl = {0,0,0};
        return InterpolatedDispl;
    }

    // Also barycentric coordinates need to be positive numbers (else the coordinate is outside of the cell).
    // So if one of the coordinates is smaller than zero
    float eps = 0.0;
//    if(BaryCoord[0] <= 0.0 || BaryCoord[1] <= 0.0 || BaryCoord[2] <= 0.0 || BaryCoord[3] <= 0.0)
//    if (BaryCoord[0] < 0.0 || BaryCoord[0] > 1.0 || BaryCoord[1] < 0.0 || BaryCoord[1] > 1.0 ||
//        BaryCoord[2] < 0.0 || BaryCoord[2] > 1.0 || BaryCoord[3] < 0.0 || BaryCoord[3] > 1.0 ) {
    if (BaryCoord[0] < 0.0 - eps || BaryCoord[0] > 1.0 + eps || BaryCoord[1] < 0.0 - eps || BaryCoord[1] > 1.0 + eps ||
        BaryCoord[2] < 0.0 - eps || BaryCoord[2] > 1.0 + eps || BaryCoord[3] < 0.0 - eps || BaryCoord[3] > 1.0 + eps) {
//    if (BaryCoord[0] <= 0.0 || BaryCoord[1] <= 0.0 || BaryCoord[2] <= 0.0 || BaryCoord[3] <= 0.0) {
        // Set displacement zero and end function immediately!
//        std::cout<<"There is negative barycentric coordinate at this D grid point! "<<std::endl;
        InterpolatedDispl = {float_max, float_max, float_max};
//        if (Map) { InterpolatedDispl = {float_max, float_max, float_max}; }
//        else { InterpolatedDispl = {0, 0, 0}; }
//        InterpolatedDispl = {0,0,0};
        return InterpolatedDispl;


}
 */