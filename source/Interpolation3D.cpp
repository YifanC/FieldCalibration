#include "../include/Interpolation3D.hpp"
#include "../include/Laser.hpp"
#include "../include/ThreeVector.hpp"
#include <limits>

std::vector<std::pair<unsigned, float>>
GetClosestTracksInfo(std::vector<LaserTrack> &LaserTrackSet, const unsigned NumberOfClosestTracks) {
//   std::vector<LaserTrack> ClosestTracks;
//   ClosestTracks.resize(NumberOfClosestTracks);

    std::vector<std::pair<unsigned, float>> ClosestTracksInfo;
    for (unsigned NumberOfTracks = 0; NumberOfTracks < NumberOfClosestTracks; NumberOfTracks++) {
        ClosestTracksInfo.push_back(std::make_pair(0, 0xDEADBEEF));
    }


    ThreeVector<float> PoyntingVector = {1.0, 1.0, 3.0};
    PoyntingVector -= LaserTrackSet.front().GetLaserPosition();
    PoyntingVector /= PoyntingVector.GetNorm();
    std::array<float, 2> AnglesOfPoint = AnglesFromPoynting(PoyntingVector);

    std::cout << AnglesOfPoint[0] * 180 / M_PI << " " << AnglesOfPoint[1] * 180 / M_PI << std::endl;

    for (unsigned track = 0; track < LaserTrackSet.size(); track++) {
        float Radius = 0;
        for (unsigned angle = 0; angle < AnglesOfPoint.size(); angle++)
            Radius += std::pow(AnglesOfPoint[angle] - LaserTrackSet[track].GetAngles()[angle], 2);

        if (ClosestTracksInfo.back().second > Radius) {
            ClosestTracksInfo.back().first = track;
            ClosestTracksInfo.back().second = Radius;

            std::sort(ClosestTracksInfo.begin(), ClosestTracksInfo.end(), PairSortFunction);
//       std::cout << ClosestTracksInfo[0].first << " " << ClosestTracksInfo[0].second << " " << ClosestTracksInfo[1].first  << " " << ClosestTracksInfo[1].second << " " 
// 		<< ClosestTracksInfo[2].first << " " << ClosestTracksInfo[2].second /*<< " " << ClosestTracksInfo[3].first  << " " << ClosestTracksInfo[3].second */<< std::endl; 
        }
    }

    return ClosestTracksInfo;
}

std::vector<std::pair<unsigned int, unsigned int>>
GetClosestLaserSample(std::vector<LaserTrack> &LaserTrackSet, const unsigned int NumberOfClosestSamples) {
    ThreeVector<float> InterpolPoint = {1.0, 1.0, 3.0};

    std::vector<std::pair<unsigned int, unsigned int>> ClosestLaserSample;
    ClosestLaserSample.resize(NumberOfClosestSamples);

    std::vector<std::pair<unsigned int, float>> ClosestSampleInfo;
    for (unsigned info = 0; info < NumberOfClosestSamples; info++)
        ClosestSampleInfo.push_back(std::make_pair(0, 0xDEADBEEF));

    std::vector<std::pair<unsigned int, float>> ClosestTracksInfo = GetClosestTracksInfo(LaserTrackSet, 3);

    float Radius = 0;
    for (unsigned info = 0; info < ClosestTracksInfo.size(); info++) {
        for (unsigned tracksample = 0;
             tracksample < LaserTrackSet[ClosestTracksInfo[info].first].GetNumberOfSamples(); tracksample++) {
            Radius = (InterpolPoint -
                      LaserTrackSet[ClosestTracksInfo[info].first].GetSamplePosition(tracksample)).GetNorm();
            if (!info && ClosestSampleInfo[1].second > Radius) {
                ClosestSampleInfo[1].first = tracksample;
                ClosestSampleInfo[1].second = Radius;

                std::sort(ClosestSampleInfo.begin(), ClosestSampleInfo.begin() + 1, PairSortFunction);
// 	std::cout << ClosestSampleInfo.front().first << " " << ClosestSampleInfo.front().second << " " << ClosestSampleInfo[1].first <<" " << ClosestSampleInfo[1].second << std::endl;
            } else if (info && ClosestSampleInfo.back().second > Radius) {
                ClosestSampleInfo.back().first = tracksample;
                ClosestSampleInfo.back().second = Radius;

                std::sort(ClosestSampleInfo.end() - 1, ClosestSampleInfo.end(), PairSortFunction);
            }
        }
    }
    return ClosestLaserSample;
}

std::array<float, 2> AnglesFromPoynting(ThreeVector<float> &Poynting) {
    Poynting /= Poynting.GetNorm();

    std::array<float, 2> ang_res;
    ang_res[0] = std::asin(Poynting[1]);
    ang_res[1] = std::asin(Poynting[0] / std::cos(ang_res[0]));
    return ang_res;
}

ThreeVector<float> PoyntingFromAngles(const std::array<float, 2> &Angles) {
    // Create unit vector of Poynting vector from spherical coordinates
    ThreeVector<float> vec_res;
    vec_res[0] = (float) (std::cos(Angles[0]) * std::sin(Angles[1])); // cos(theta)*cos(phi)
    vec_res[1] = (float) std::sin(Angles[0]); // sin(theta)
    vec_res[2] = (float) (std::cos(Angles[0]) * std::cos(Angles[1])); // cos(theta)*cos(phi)

    return vec_res;
}

bool PairSortFunction(std::pair<unsigned, float> left_pair, std::pair<unsigned, float> right_pair) {
    return (left_pair.second < right_pair.second);
}

// This produces the Delaunay Mesh from LaserTrackSet in a std::vector
Delaunay TrackMesher(const std::vector<LaserTrack> &LaserTrackSet) {
    // Create a vector with a std::pair< Point, PointInfo (also a pair) > this is the input format for the Delaunay constructor
    std::vector<std::pair<Point, std::pair<unsigned long, unsigned long>>> Points;

    // Loop over tracks in the vector
    for (unsigned long track = 0; track < LaserTrackSet.size(); track++) {
        // Loop over data points (samples) of each track
        for (unsigned long sample = 0; sample < LaserTrackSet[track].GetNumberOfSamples(); sample++) {
            // Convert ThreeVector<float> of sample position to CGAL point
            Point SamplePoint = Point(LaserTrackSet[track].GetSamplePosition(sample)[0],
                                      LaserTrackSet[track].GetSamplePosition(sample)[1],
                                      LaserTrackSet[track].GetSamplePosition(sample)[2]);

            // Prepare point info (pair with track number and track sample number)
            std::pair<unsigned long, unsigned long> SamplePointIndex = std::make_pair(track, sample);

            // Fill Delaunay input container with point and point info
            Points.push_back(std::make_pair(SamplePoint, SamplePointIndex));
        } // end sample loop
    } // end track loop

    // Create Delaunay mesh (this takes quite some runtime!)
    Delaunay DelaunayMesh(Points.begin(), Points.end());

    // Return mesh
    return DelaunayMesh;
}

//Attention! the typedef of EMAP related Delaunay template is with x in front
xDelaunay Mesher(std::vector<ThreeVector<float>> &Position, TPCVolumeHandler &TPC) {
    // Create a vector with a std::pair< Point, PointInfo (also a pair) > this is the input format for the Delaunay constructor
    std::vector<std::pair<Point, int> > Points;

    // Loop over the size of grid
    for (int n = 0; n < Position.size(); n++) {
        // Convert ThreeVector<float> of sample position to CGAL point
        xPoint SamplePoint = xPoint(Position[n][0], Position[n][1], Position[n][2]);

        // Prepare point info is index of the mesh points
        int SamplePointIndex = n;

        // Fill Delaunay input container with point and point info
        Points.push_back(std::make_pair(SamplePoint, SamplePointIndex));

    }

    // Create Delaunay mesh (this takes quite some runtime!)
    xDelaunay DelaunayMesh(Points.begin(), Points.end());

    // Return mesh
    return DelaunayMesh;
}

ThreeVector<float> PointToVector(Point &InputPoint) {
    ThreeVector<float> vec_res = {(float) InputPoint[0], (float) InputPoint[1], (float) InputPoint[2]};
    return vec_res;
}

Point VectorToPoint(ThreeVector<float> &InputVector) {
    Point point_res(InputVector[0], InputVector[1], InputVector[2]);
    return point_res;
}

//////////Is this necessary?
xPoint xVectorToPoint(ThreeVector<float> &InputVector) {
    xPoint point_res(InputVector[0], InputVector[1], InputVector[2]);
    return point_res;
}

// This function Interpolates the displacement of Location within the Mesh
// For InterpolateCGAL/InterpolateMap,
// the first parameter is the one to provide the vector (e.g displacement)
// the second parameter is the one to indicate the mesh position
// If the Map flag is true: set the unknown point to float max, which is suitable for the final map establish
// If the Map flag is false: set the unknown point to 0, which is suitable for the middle steps
ThreeVector<float>
InterpolateCGAL(const std::vector<LaserTrack> &LaserTrackSet, const std::vector<LaserTrack> &LaserMeshSet,
                const Delaunay &Mesh, ThreeVector<float> Location) {
    float float_max = std::numeric_limits<float>::max();

    // Create a array which contains the info of all 4 vertices of a cell
    std::array<std::pair<unsigned long, unsigned long>, 4> PointIndex;

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

    // If the function is still alive, loop over all barycentric coordinates
    for (unsigned vertex_no = 0; vertex_no < 4; vertex_no++) {
        // Use the barycentric coordinates as a weight for the correction stored at this vertex in order to get the interpolated displacement
        InterpolatedDispl += (LaserTrackSet[PointIndex[vertex_no].first].GetDisplacement(PointIndex[vertex_no].second) *
                              BaryCoord[vertex_no]);
    }

    // Return interpolated displacement
    return InterpolatedDispl;
}

////Position refers to the true points which are deduced from the reconstructed grid of Correction Map
////Mesh refers to the result of Mesher()
////Location is the new grid point(true space coordinate)
ThreeVector<float> EInterpolateCGAL(std::vector<ThreeVector<float>> &En, std::vector<ThreeVector<float>> &Position,
                                    const xDelaunay &Mesh, ThreeVector<float> Location, const TPCVolumeHandler &TPC) {
    float float_max = std::numeric_limits<float>::max();

    ThreeVector<unsigned long> Reso = TPC.GetDetectorResolution();

    // Create a array which contains the info of all 4 vertices of a cell
    std::array<int, 4> Index;

    // Initialize a displacement vector with zero
    ThreeVector<float> InterpolatedEfield = {0.0, 0.0, 0.0};

    // Initialize Barycentric coordinate system (it will have 4 dimensions)
    std::vector<float> BaryCoord;

    // Find cell in the mesh where the point is located
//    xDelaunay::Cell_handle Cell =  Mesh.locate(xVectorToPoint(Location));

    xDelaunay::Locate_type loc;
    int li, lj;
    xDelaunay::Cell_handle Cell = Mesh.locate(xVectorToPoint(Location), loc, li, lj);
//    std::cout<<"loc: "<<loc<<"; li: "<<li<<"; lj: "<<lj<<std::endl;

    // Loop over all four vertex points of the cell of interest
    for (unsigned vertex_no = 0; vertex_no < Index.size(); vertex_no++) {
        // Get vertex info "n" (the index of the Position vector [Attention! It must be corresponding to the index of En vector!])
        Index[vertex_no] = Cell->vertex(vertex_no)->info();
    }

    // Initialize matrix for Location transformation into barycentric coordinate system
    Matrix3x3 TransMatrix = {{0, 0, 0},
                             {0, 0, 0},
                             {0, 0, 0}};

    // Loop over transverse matrix rows and columns
    for (unsigned row = 0; row < 3; row++) //x,y,z
    {
        for (unsigned column = 0; column < 3; column++)//1,2,3
        {
            // Fill transformation matrix elements
            // Valid for 3D barycentric coordinate system
            // When loop the E local position the order is z y x, be careful of the index of the vector
            TransMatrix[row][column] = Position[Index[column]][row] - Position[Index[3]][row];
        }
    }

    // r - r4 in threevector
//    Location -= Position[Index[3]];// Reuse Location and store its position relative to the last vertex of the cell it is contained in
    ThreeVector<float> RR4 = Location - Position[Index[3]];

    // If the transformation matrix can be successfully inverted
    if (TransMatrix.Invert()) {
        // Use inverted matrix to fill the first three coordinates
        ThreeVector<float> BC = TransMatrix * RR4;
        BaryCoord = BC.GetStdVector();

        // The sum of all barycentric coordinates has to be 1 by definition, use this to calculate the 4th coordinate
        BaryCoord.push_back(1 - BaryCoord[0] - BaryCoord[1] - BaryCoord[2]);
    } else // if the matrix can't be inverted
    {
        // Set E field to zero and end function immediately!
//        std::cout<<"The transition matrix for this E grid point is not invertable. "<<std::endl;
//        InterpolatedEfield = {273.0,0.0,0.0};
        InterpolatedEfield = {float_max, float_max, float_max};
//        InterpolatedEfield = {-99,-99,-99};
        return InterpolatedEfield;
    }

    // Also barycentric coordinates need to be positive numbers (else the coordinate is outside of the cell).
    // So if one of the coordinates is negative, terminate the function
//    float eps = 1E4;
    float eps = 0.0;
//    if(BaryCoord[0] <= 0.0 || BaryCoord[1] <= 0.0 || BaryCoord[2] <= 0.0 || BaryCoord[3] <= 0.0)
    if (BaryCoord[0] < 0.0 - eps || BaryCoord[0] > 1.0 + eps || BaryCoord[1] < 0.0 - eps || BaryCoord[1] > 1.0 + eps ||
        BaryCoord[2] < 0.0 - eps || BaryCoord[2] > 1.0 + eps || BaryCoord[3] < 0.0 - eps || BaryCoord[3] > 1.0 + eps) {
        // Set E field to zero and end function immediately!
//        std::cout<<"There is negative barycentric coordinate at this E grid point! "<<std::endl;
//        InterpolatedEfield = {273.0,0.0,0.0};
//        std::cout << "x: " << Location[0]<<"; y: "<<Location[1]<<"; z: "<<Location[2]<<std::endl;
//        std::cout<<"A: "<< BaryCoord[0]<<"; B: "<< BaryCoord[1]<<"; C: "<< BaryCoord[2]<<"; D: "<< BaryCoord[3]<<std::endl;
//        InterpolatedEfield = {-99,-99,-99};
//            std::cout<<"loc: "<<loc<<"; li: "<<li<<"; lj: "<<lj<<std::endl;
        InterpolatedEfield = {float_max, float_max, float_max};
        return InterpolatedEfield;
    }

    // If the function is still alive, loop over all barycentric coordinates
    for (unsigned vertex_no = 0; vertex_no < 4; vertex_no++) {
        // Use the barycentric coordinates as a weight for the correction stored at this vertex in order to get the interpolated displacement
        // Adding up the barycoord components as a whole vector of (x1,y1,z1) instead of x1,x2,x3....
        // BaryCoord[vertex_no] is a number
        InterpolatedEfield += En[Index[vertex_no]] * BaryCoord[vertex_no];
    }

    // Return interpolated E field
    if(InterpolatedEfield[0]>20){
        std::cout<<"Ex: "<< InterpolatedEfield[0] <<"; Ey: "<<InterpolatedEfield[1]<<"; Ez: "<<InterpolatedEfield[2]<<std::endl;
    }
    return InterpolatedEfield;
}

//////Position refers to the true points which are deduced from the reconstructed grid of Correction Map
//////Mesh refers to the result of Mesher()
//////Location is the new grid point(true space coordinate)
//float EcompInterpolateCGAL(std::vector<float> &En, std::vector<ThreeVector<float>> &Position,
//                                        const xDelaunay &Mesh, ThreeVector<float> Location, const TPCVolumeHandler &TPC) {
//    float float_max = std::numeric_limits<float>::max();
//
//    ThreeVector<unsigned long> Reso = TPC.GetDetectorResolution();
//
//    // Create a array which contains the info of all 4 vertices of a cell
//    std::array<int, 4> Index;
//
//    // Initialize a displacement vector with zero
//    ThreeVector<float> InterpolatedEfield = {0.0, 0.0, 0.0};
//
//    // Initialize Barycentric coordinate system (it will have 4 dimensions)
//    std::vector<float> BaryCoord;
//
//    // Find cell in the mesh where the point is located
////    xDelaunay::Cell_handle Cell =  Mesh.locate(xVectorToPoint(Location));
//
//    xDelaunay::Locate_type loc;
//    int li, lj;
//    xDelaunay::Cell_handle Cell = Mesh.locate(xVectorToPoint(Location), loc, li, lj);
////    std::cout<<"loc: "<<loc<<"; li: "<<li<<"; lj: "<<lj<<std::endl;
//
//    // Loop over all four vertex points of the cell of interest
//    for (unsigned vertex_no = 0; vertex_no < Index.size(); vertex_no++) {
//        // Get vertex info "n" (the index of the Position vector [Attention! It must be corresponding to the index of En vector!])
//        Index[vertex_no] = Cell->vertex(vertex_no)->info();
//    }
//
//    // Initialize matrix for Location transformation into barycentric coordinate system
//    Matrix3x3 TransMatrix = {{0, 0, 0},
//                             {0, 0, 0},
//                             {0, 0, 0}};
//
//    // Loop over transverse matrix rows and columns
//    for (unsigned row = 0; row < 3; row++) //x,y,z
//    {
//        for (unsigned column = 0; column < 3; column++)//1,2,3
//        {
//            // Fill transformation matrix elements
//            // Valid for 3D barycentric coordinate system
//            // When loop the E local position the order is z y x, be careful of the index of the vector
//            TransMatrix[row][column] = Position[Index[column]][row] - Position[Index[3]][row];
//        }
//    }
//
//    // r - r4 in threevector
////    Location -= Position[Index[3]];// Reuse Location and store its position relative to the last vertex of the cell it is contained in
//    ThreeVector<float> RR4 = Location - Position[Index[3]];
//
//    // If the transformation matrix can be successfully inverted
//    if (TransMatrix.Invert()) {
//        // Use inverted matrix to fill the first three coordinates
//        ThreeVector<float> BC = TransMatrix * RR4;
//        BaryCoord = BC.GetStdVector();
//
//        // The sum of all barycentric coordinates has to be 1 by definition, use this to calculate the 4th coordinate
//        BaryCoord.push_back(1 - BaryCoord[0] - BaryCoord[1] - BaryCoord[2]);
//    } else // if the matrix can't be inverted
//    {
//        // Set E field to zero and end function immediately!
////        std::cout<<"The transition matrix for this E grid point is not invertable. "<<std::endl;
////        InterpolatedEfield = {273.0,0.0,0.0};
//        InterpolatedEfield = {float_max, float_max, float_max};
////        InterpolatedEfield = {-99,-99,-99};
//        return InterpolatedEfield;
//    }
//
//    // Also barycentric coordinates need to be positive numbers (else the coordinate is outside of the cell).
//    // So if one of the coordinates is negative, terminate the function
////    float eps = 1E4;
//    float eps = 0.0;
////    if(BaryCoord[0] <= 0.0 || BaryCoord[1] <= 0.0 || BaryCoord[2] <= 0.0 || BaryCoord[3] <= 0.0)
//    if (BaryCoord[0] < 0.0 - eps || BaryCoord[0] > 1.0 + eps || BaryCoord[1] < 0.0 - eps || BaryCoord[1] > 1.0 + eps ||
//        BaryCoord[2] < 0.0 - eps || BaryCoord[2] > 1.0 + eps || BaryCoord[3] < 0.0 - eps || BaryCoord[3] > 1.0 + eps) {
//        // Set E field to zero and end function immediately!
////        std::cout<<"There is negative barycentric coordinate at this E grid point! "<<std::endl;
////        InterpolatedEfield = {273.0,0.0,0.0};
////        std::cout << "x: " << Location[0]<<"; y: "<<Location[1]<<"; z: "<<Location[2]<<std::endl;
////        std::cout<<"A: "<< BaryCoord[0]<<"; B: "<< BaryCoord[1]<<"; C: "<< BaryCoord[2]<<"; D: "<< BaryCoord[3]<<std::endl;
////        InterpolatedEfield = {-99,-99,-99};
////            std::cout<<"loc: "<<loc<<"; li: "<<li<<"; lj: "<<lj<<std::endl;
//        InterpolatedEfield = {float_max, float_max, float_max};
//        return InterpolatedEfield;
//    }
//
//    // If the function is still alive, loop over all barycentric coordinates
//    for (unsigned vertex_no = 0; vertex_no < 4; vertex_no++) {
//        // Use the barycentric coordinates as a weight for the correction stored at this vertex in order to get the interpolated displacement
//        // Adding up the barycoord components as a whole vector of (x1,y1,z1) instead of x1,x2,x3....
//        // BaryCoord[vertex_no] is a number
//        InterpolatedEfield += En[Index[vertex_no]] * BaryCoord[vertex_no];
//    }
//
//    // Return interpolated E field
//    return InterpolatedEfield;
//}


// This function interpolates regularly spaced grid points of the TPC and stores them in a std::vector (can later be used in the WriteRootFile function)
// For InterpolateCGAL/InterpolateMap,
// the first parameter is the one to provide the vector (e.g displacement)
// the second parameter is the one to indicate the mesh position
// Set CorrMapFlag to false as default in head file
std::vector<ThreeVector<float>>
InterpolateMap(const std::vector<LaserTrack> &LaserTrackSet, const std::vector<LaserTrack> &LaserMeshSet,
               const Delaunay &Mesh, const TPCVolumeHandler &TPC) {
    // Initialize output data structure
    std::vector<ThreeVector<float>> DisplacementMap;

    // Initialize temporary location vector
    ThreeVector<float> Location;
    
    // Define the bin size of TPC volume
    ThreeVector<float> Unit = {TPC.GetDetectorSize()[0] / static_cast<float>(TPC.GetDetectorResolution()[0] - 1),
                               TPC.GetDetectorSize()[1] / static_cast<float>(TPC.GetDetectorResolution()[1] - 1),
                               TPC.GetDetectorSize()[2] / static_cast<float>(TPC.GetDetectorResolution()[2] - 1)};

    // Loop over all xbins of the TPC
    for (unsigned xbin = 0; xbin < TPC.GetDetectorResolution()[0]; xbin++) {
        std::cout << "Processing plane " << xbin << " of " << TPC.GetDetectorResolution()[0] << std::endl;
        // Calculate Grid point x-coordinate
        Location[0] = TPC.GetDetectorOffset()[0] + Unit[0] * xbin;

        // Loop over all ybins of the TPC
        for (unsigned ybin = 0; ybin < TPC.GetDetectorResolution()[1]; ybin++) {
            // Calculate Grid point y-coordinate
            Location[1] = TPC.GetDetectorOffset()[1] + Unit[1] * ybin;

            // Loop over all zbins of the TPC
            for (unsigned zbin = 0; zbin < TPC.GetDetectorResolution()[2]; zbin++) {
                // Calculate Grid point y-coordinate
                Location[2] = TPC.GetDetectorOffset()[2] + Unit[2] * zbin;

                // Fill displacement map
                // The Map trigger is turned on for interpolation for regular grid. This will set the unknown point to float_max
                DisplacementMap.push_back(InterpolateCGAL(LaserTrackSet, LaserMeshSet, Mesh, Location));

            } // end zbin loop
        } // end ybin loop
    } // end ybin loop

    return DisplacementMap;
}

//// This part can be reduced with template maybe..
//// This function interpolates regularly spaced grid points of the TPC and stores them in a std::vector (can later be used in the WriteRootFile function)
std::vector<ThreeVector<float>>
EInterpolateMap(std::vector<ThreeVector<float>> &En, std::vector<ThreeVector<float>> &Position, const xDelaunay &Mesh,
                const TPCVolumeHandler &TPC, ThreeVector<unsigned long> EReso) {
    // Initialize output data structure
    std::vector<ThreeVector<float>> EMap;

    // Initialize temporary location vector
    ThreeVector<float> Location;

    // Loop over all xbins of the TPC
//    for(unsigned xbin = 0; xbin < TPC.GetDetectorResolution()[0]; xbin++)
    for (unsigned xbin = 0; xbin < EReso[0]; xbin++) {
        std::cout << "Processing plane " << xbin + 1 << " of " << EReso[0] << " for Emap. " << std::endl;

        // Calculate Grid point x-coordinate
        Location[0] = TPC.GetDetectorOffset()[0] + TPC.GetDetectorSize()[0] / static_cast<float>(EReso[0] - 1) * xbin;

        // Loop over all ybins of the TPC
        for (unsigned ybin = 0; ybin < EReso[1]; ybin++) {
            // Calculate Grid point y-coordinate
            Location[1] =
                    TPC.GetDetectorOffset()[1] + TPC.GetDetectorSize()[1] / static_cast<float>(EReso[1] - 1) * ybin;

            // Loop over all zbins of the TPC
            for (unsigned zbin = 0; zbin < EReso[2]; zbin++) {
                // Calculate Grid point y-coordinate
                Location[2] =
                        TPC.GetDetectorOffset()[2] + TPC.GetDetectorSize()[2] / static_cast<float>(EReso[2] - 1) * zbin;

//                std::cout<<"x: "<<Location[0]<<"; y: "<<Location[1]<<"; z: "<<Location[2]<<std::endl;

                // Fill displacement map
                EMap.push_back(EInterpolateCGAL(En, Position, Mesh, Location, TPC));
            } // end zbin loop
        } // end ybin loop
    } // end xbin loop

    return EMap;
}


void InterpolateTrack(LaserTrack &Track, const std::vector<LaserTrack> &LaserTrackSet, const Delaunay &Mesh) {
    std::array<std::pair<unsigned long, unsigned long>, 4> PointIndex;
    ThreeVector<float> InterpolatedDispl = {0.0, 0.0, 0.0};
    ThreeVector<float> Location;

    for (unsigned sample_no = 0; sample_no < Track.GetNumberOfSamples(); sample_no++) {
        Location = Track.GetSamplePosition(sample_no);
        Delaunay::Cell_handle Cell = Mesh.locate(VectorToPoint(Location));

        ThreeVector<float> VertexPoint;

        for (unsigned vertex_no = 0; vertex_no < 4; vertex_no++) {
            PointIndex[vertex_no] = Cell->vertex(vertex_no)->info();
        }

        Matrix3x3 TransMatrix = {{0, 0, 0},
                                 {0, 0, 0},
                                 {0, 0, 0}};

        for (unsigned row = 0; row < 3; row++) {
            for (unsigned column = 0; column < 3; column++) {
                TransMatrix[row][column] =
                        LaserTrackSet[PointIndex[column].first].GetSamplePosition(PointIndex[column].second)[row] -
                        LaserTrackSet[PointIndex.back().first].GetSamplePosition(PointIndex.back().second)[row];
            }
        }

        Location -= LaserTrackSet[PointIndex.back().first].GetSamplePosition(PointIndex.back().second);
//     std::cout << TransMatrix.Determinant() << std::endl;

        std::vector<float> BaryCoord;
        if (TransMatrix.Invert()) {
            BaryCoord = (TransMatrix * Location).GetStdVector();
            BaryCoord.push_back(1 - BaryCoord[0] - BaryCoord[1] - BaryCoord[2]);
//     std::cout << BaryCoord[0] << " " << BaryCoord[1] << " " << BaryCoord[2] << " " << BaryCoord[3] << std::endl;
        } else {
//       BaryCoord = {0.0,0.0,0.0,0.0};
            InterpolatedDispl = {0.0, 0.0, 0.0};
            return;
        }

        if (BaryCoord[0] < 0.0 || BaryCoord[1] < 0.0 || BaryCoord[2] < 0.0 || BaryCoord[3] < 0.0) {
//       std::cout << BaryCoord[0] << " " << BaryCoord[1] << " " << BaryCoord[2] << " " << BaryCoord[3] << std::endl;
            InterpolatedDispl = {0.0, 0.0, 0.0};
            return;
        }

        for (unsigned vertex_no = 0; vertex_no < 4; vertex_no++) {
            InterpolatedDispl += (
                    LaserTrackSet[PointIndex[vertex_no].first].GetDisplacement(PointIndex[vertex_no].second) *
                    BaryCoord[vertex_no]);
        }

        InterpolatedDispl = ThreeVector<float>::DotProduct(InterpolatedDispl, Track.GetPoyntingVector()) *
                            Track.GetPoyntingVector();
        Track.AddToDisplacement(InterpolatedDispl, sample_no);
    }
}
