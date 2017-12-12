//
// Created by matthias on 28.04.17.
//

#include <TVector3.h>

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../include/Utilities.hpp"
#include "../include/LaserTrack.hpp"
#include "../include/Laser.hpp"
#include "../include/Matrix3x3.hpp"
#include "../include/ThreeVector.hpp"
#include "../include/Interpolation3D.hpp"

// Tests factorial of 0.
TEST(TestSplitter, CheckSelection) {

    TVector3 entry(1.,1.,1.);
    TVector3 exit(0.,0.,0.);

    std::vector<TVector3> RecobLaserTrack;

    // fill the dummy vector in the right order. Tracks are assumed to be filled ordered,
    // starting at the entry point and ending closest to the exit point.
    for (int i=9; i > 0; i--){
        const float pt = i / 10.0;
        RecobLaserTrack.push_back(TVector3(pt,pt,pt));
    }

    // Make a laser track out of entry, exit and track points. Then add it to the Laser collection
    LaserTrack Track1 = LaserTrack(entry, exit, RecobLaserTrack);
    std::vector<LaserTrack> PewPew{ Track1 };
    Laser Las(PewPew);

    // Now the actual testing of the function follows
    // Check if track gets into the first selection:
    auto res = ReachedExitPoint(PewPew, 0.5);
    auto InSelection = res.front();
    auto OutSelection = res.back();
    ASSERT_EQ(InSelection.GetNumberOfTracks(), 1);
    ASSERT_EQ(OutSelection.GetNumberOfTracks(), 0);
    res.clear();

    // Check if track gets into the second selection under new condition:
    res = ReachedExitPoint(PewPew, 0.001);
    InSelection = res.front();
    OutSelection = res.back();
    ASSERT_EQ(InSelection.GetNumberOfTracks(), 0);
    ASSERT_EQ(OutSelection.GetNumberOfTracks(), 1);

    //float a = 1;
    ASSERT_TRUE(true);
    //EXPECT_EQ(1, ReachedExitPoint(Test, 1));
}

TEST(TestDownsampler, SplitSingleEvenTrack) {

    TVector3 entry(1.,1.,1.);
    TVector3 exit(0.,0.,0.);

    std::vector<TVector3> RecobLaserTrack;

    // fill the dummy vector in the right order. Tracks are assumed to be filled ordered,
    // starting at the entry point and ending closest to the exit point.
    for (int i=0; i < 10; i++){
        const float pt = i;
        RecobLaserTrack.push_back(TVector3(pt,pt,pt));
    }

    // Make a laser track out of entry, exit and track points. Then add it to the Laser collection
    LaserTrack Track1 = LaserTrack(entry, exit, RecobLaserTrack);
    std::vector<LaserTrack> PewPew{ Track1 };
    Laser Las(PewPew);

    std::vector<Laser> LaserSets = SplitTrackSet(Las, 2);

    auto first_set = LaserSets.front();
    auto first_track = first_set.GetFirstTrack();

    ASSERT_TRUE(first_track.GetReco().front() == ThreeVector<float>(0., 0., 0.));
    ASSERT_TRUE(first_track.GetReco().back() == ThreeVector<float>(8., 8., 8.));
    ASSERT_EQ(first_track.GetNumberOfSamples(), 5);

    auto second_set = LaserSets.back();
    auto second_track = second_set.GetFirstTrack();

    ASSERT_TRUE(second_track.GetReco().front() == ThreeVector<float>(1., 1., 1.));
    ASSERT_TRUE(second_track.GetReco().back() == ThreeVector<float>(9., 9., 9.));
    ASSERT_EQ(second_track.GetNumberOfSamples(), 5);
}

TEST(TestDownsampler, SplitSingleOddTrack) {

    TVector3 entry(1.,1.,1.);
    TVector3 exit(0.,0.,0.);

    std::vector<TVector3> RecobLaserTrack;

    // fill the dummy vector in the right order. Tracks are assumed to be filled ordered,
    // starting at the entry point and ending closest to the exit point.
    for (int i=0; i < 11; i++){
        const float pt = i;
        RecobLaserTrack.push_back(TVector3(pt,pt,pt));
    }

    // Make a laser track out of entry, exit and track points. Then add it to the Laser collection
    LaserTrack Track1 = LaserTrack(entry, exit, RecobLaserTrack);
    std::vector<LaserTrack> PewPew{ Track1 };
    Laser Las(PewPew);

    std::vector<Laser> LaserSets = SplitTrackSet(Las, 2);

    auto first_set = LaserSets.front();
    auto first_track = first_set.GetFirstTrack();

    ASSERT_TRUE(first_track.GetReco().front() == ThreeVector<float>(0., 0., 0.));
    ASSERT_TRUE(first_track.GetReco().back() == ThreeVector<float>(10., 10., 10.));
    ASSERT_EQ(first_track.GetNumberOfSamples(), 6);

    auto second_set = LaserSets.back();
    auto second_track = second_set.GetFirstTrack();

    ASSERT_TRUE(second_track.GetReco().front() == ThreeVector<float>(1., 1., 1.));
    ASSERT_TRUE(second_track.GetReco().back() == ThreeVector<float>(9., 9., 9.));
    ASSERT_EQ(second_track.GetNumberOfSamples(), 5);
}

TEST(TestDownsampler, HigherOrder) {

    TVector3 entry(1.,1.,1.);
    TVector3 exit(0.,0.,0.);

    std::vector<TVector3> RecobLaserTrack;

    // fill the dummy vector in the right order. Tracks are assumed to be filled ordered,
    // starting at the entry point and ending closest to the exit point.
    for (int i=1; i < 11; i++){
        const float pt = i;
        RecobLaserTrack.push_back(TVector3(pt,pt,pt));
    }

    // Make a laser track out of entry, exit and track points. Then add it to the Laser collection
    LaserTrack Track1 = LaserTrack(entry, exit, RecobLaserTrack);
    std::vector<LaserTrack> PewPew{ Track1 };
    Laser Las(PewPew);

    std::vector<Laser> LaserSets = SplitTrackSet(Las, 4);

    ASSERT_EQ(LaserSets.size(), 4);

    auto first_set = LaserSets[0];
    auto first_track = first_set.GetFirstTrack();

    ASSERT_TRUE(first_track.GetReco().front() == ThreeVector<float>(1., 1., 1.));
    ASSERT_TRUE(first_track.GetReco()[1] == ThreeVector<float>(5., 5., 5.));
    ASSERT_TRUE(first_track.GetReco().back() == ThreeVector<float>(9., 9., 9.));
    ASSERT_EQ(first_track.GetNumberOfSamples(), 3);

    auto second_set = LaserSets.back();
    auto second_track = second_set.GetFirstTrack();

    ASSERT_TRUE(second_track.GetReco().front() == ThreeVector<float>(4., 4., 4.));
    ASSERT_TRUE(second_track.GetReco().back() == ThreeVector<float>(8., 8., 8.));
    ASSERT_EQ(second_track.GetNumberOfSamples(), 2);

    ASSERT_EQ(LaserSets[1].GetFirstTrack().GetNumberOfSamples(), 3);
    ASSERT_EQ(LaserSets[2].GetFirstTrack().GetNumberOfSamples(), 2);
}

//TEST(TestDownsampler, TwoTracks) {
//
//    TVector3 entry(1.,1.,1.);
//    TVector3 exit(0.,0.,0.);
//
//    std::vector<TVector3> RecobLaserTrack1;
//    std::vector<TVector3> RecobLaserTrack2;
//
//    // fill the dummy vector in the right order. Tracks are assumed to be filled ordered,
//    // starting at the entry point and ending closest to the exit point.
//    for (int i=1; i < 11; i++){
//        const float pt = i;
//        RecobLaserTrack1.push_back(TVector3(pt,pt,pt));
//        RecobLaserTrack2.push_back(TVector3(-pt,-pt,-pt));
//    }
//
//    // Make a laser track out of entry, exit and track points. Then add it to the Laser collection
//    LaserTrack Track1 = LaserTrack(entry, exit, RecobLaserTrack1);
//    LaserTrack Track2 = LaserTrack(entry, exit, RecobLaserTrack2);
//    std::vector<LaserTrack> PewPew{ Track1, Track2 };
//    Laser Las(PewPew);
//
//    std::vector<Laser> LaserSets = SplitTrackSet(Las, 2);
//
//    ASSERT_EQ(LaserSets.size(), 2);
//
//    auto first_set = LaserSets[0];
//    auto second_set = LaserSets[1];
//
//    ASSERT_EQ( first_set.GetNumberOfTracks(), 2);
//    ASSERT_EQ(second_set.GetNumberOfTracks(), 2);
//
//    auto reco_set0_track0 = first_set.GetTrack(0);
//    auto reco_set0_track1 = first_set.GetTrack(1);
//
//    auto reco_set1_track0 = second_set.GetTrack(0);
//    auto reco_set1_track1 = second_set.GetTrack(1);
//
//
//    ASSERT_TRUE(reco_set0_track0.GetReco().front() == ThreeVector<float>(1., 1., 1.));
//    ASSERT_TRUE(reco_set1_track0.GetReco().front() == ThreeVector<float>(2., 2., 2.));
//
//    ASSERT_TRUE(reco_set0_track1.GetReco().front() == ThreeVector<float>(-1., -1., -1.));
//    ASSERT_TRUE(reco_set1_track1.GetReco().front() == ThreeVector<float>(-2., -2., -2.));
//
//}

TEST(Interpolation, BaryCentric) {

    float float_max = std::numeric_limits<float>::max();

    std::vector<ThreeVector<float>> Vertex(4, ThreeVector<float>(0., 0., 0.));

    Vertex[0] = {0.0, 0.0, 0.0};
    Vertex[1] = {1.0, 0.0, 0.0};
    Vertex[2] = {0.0, 1.0, 0.0};
    Vertex[3] = {0.0, 0.0, 1.0};

    std::vector<ThreeVector<float>> DisplVector(4, ThreeVector<float>(0., 0., 0.));

    DisplVector[0] = {1.0, 2.0, 3.0};
    DisplVector[1] = {4.0, 5.0, 6.0};
    DisplVector[2] = {7.0, 8.0, 9.0};
    DisplVector[3] = {10.0, 11.0, 12.0};


    ThreeVector<float> Location = {0.1, 0.1, 0.1};

//    ASSERT_TRUE(Location == ThreeVector<float>(0.5, 0.5, 0.5));

// Create a array which contains the info of all 4 vertices of a cell
std::array<std::pair<unsigned long, unsigned long>, 4> PointIndex;

// Initialize a displacement vector with zero
ThreeVector<float> InterpolatedDispl = {0.0, 0.0, 0.0};

// Initialize Barycentric coordinate system (it will have 4 dimensions)
std::vector<float> BaryCoord;

//// Find cell in the mesh where the point is located
//Delaunay::Cell_handle Cell = Mesh.locate(VectorToPoint(Location));
//
//// Loop over all four vertex points of the cell of interest
//for (unsigned vertex_no = 0; vertex_no < PointIndex.size(); vertex_no++) {
//// Get vertex info of the cell (track number, sample number)
//PointIndex[vertex_no] = Cell->vertex(vertex_no)->info();
//}

    // Initialize matrix for Location transformation into barycentric coordinate system
    Matrix3x3 TransMatrix = {{0, 0, 0},
                             {0, 0, 0},
                             {0, 0, 0}};

    // Loop over matrix rows
    for (unsigned row = 0; row < 3; row++) {
        // Loop over matrix columns
        for (unsigned column = 0; column < 3; column++) {
            // Fill transformation matrix elements
            TransMatrix[row][column] = Vertex[column][row] - Vertex[3][row];
//            ASSERT_EQ(TransMatrix[row][column],)

//            LaserMeshSet[PointIndex[column].first].GetSamplePosition(PointIndex[column].second)[row] -
//            LaserMeshSet[PointIndex.back().first].GetSamplePosition(PointIndex.back().second)[row];
        }
    }

    ASSERT_EQ(TransMatrix[0][0],0.);
    ASSERT_EQ(TransMatrix[0][1],1.);
    ASSERT_EQ(TransMatrix[0][2],0.);
    ASSERT_EQ(TransMatrix[1][0],0.);
    ASSERT_EQ(TransMatrix[1][1],0.);
    ASSERT_EQ(TransMatrix[1][2],1.);
    ASSERT_EQ(TransMatrix[2][0],-1.);
    ASSERT_EQ(TransMatrix[2][1],-1.);
    ASSERT_EQ(TransMatrix[2][2],-1.);




//    ThreeVector<float> RR4 = Location - Vertex[3];

    //// Reuse Location and store its position relative to the last vertex of the cell it is contained in
    //Location -= LaserMeshSet[PointIndex.back().first].GetSamplePosition(PointIndex.back().second);
    Location -= Vertex[3];

    ASSERT_TRUE(Location == ThreeVector<float>(0.1, 0.1, -0.9));

    // If the transformation matrix can be successfully inverted
    if (TransMatrix.Invert()) {
        // Use inverted matrix to fill the first three coordinates
        ThreeVector<float> BC = TransMatrix * Location;
        BaryCoord = BC.GetStdVector();

//        ASSERT_NEAR(BC[0],1.,1E-3);
//        ASSERT_NEAR(BC[1],0.,1E-3);
//        ASSERT_NEAR(BC[2],0.,1E-3);



        // The sum of all barycentric coordinates has to be 1 by definition, use this to calculate the 4th coordinate
        BaryCoord.push_back(1 - BaryCoord[0] - BaryCoord[1] - BaryCoord[2]);

        ASSERT_EQ(BC[0],0.1);
        ASSERT_EQ(BC[1],-0.9);
        ASSERT_EQ(BC[2],0.7);
        ASSERT_EQ(BC[3],1.1);
    }
    else // if the matrix can't be inverted
    {
        //SAY SOMETHING IS WRONG!
        InterpolatedDispl = {-99.,-99.,-99.};

    }



//    ASSERT_NEAR(BaryCoord[0],1.,1E-3);
//    ASSERT_NEAR(BaryCoord[1],0.,1E-3);
//    ASSERT_NEAR(BaryCoord[2],0.,1E-3);
//    ASSERT_NEAR(BaryCoord[3],0.,1E-3);

    // Also barycentric coordinates need to be positive numbers (else the coordinate is outside of the cell).
    // So if one of the coordinates is smaller than zero
    if (BaryCoord[0] < 0.0 || BaryCoord[1] < 0.0 || BaryCoord[2] < 0.0 || BaryCoord[3] < 0.0) {
        InterpolatedDispl = {-999.,-999.,-999.};

    }

    // If the function is still alive, loop over all barycentric coordinates
    for (unsigned vertex_no = 0; vertex_no < 4; vertex_no++) {
        // Use the barycentric coordinates as a weight for the correction stored at this vertex in order to get the interpolated displacement
//        InterpolatedDispl += (LaserTrackSet[PointIndex[vertex_no].first].GetDisplacement(PointIndex[vertex_no].second) *
//        BaryCoord[vertex_no]);

        InterpolatedDispl += DisplVector[vertex_no] * BaryCoord[vertex_no];
    }

//    ASSERT_TRUE(InterpolatedDispl == ThreeVector<float>(1., 1., 1.));

//    ASSERT_NEAR(InterpolatedDispl[0],1.,1E-3);
//    ASSERT_NEAR(InterpolatedDispl[1],1.,1E-3);
//    ASSERT_NEAR(InterpolatedDispl[2],1.,1E-3);

//    ASSERT_EQ(InterpolatedDispl[0],7.);
//    ASSERT_EQ(InterpolatedDispl[1],8.);
//    ASSERT_EQ(InterpolatedDispl[2],9.);



}

TEST(Interpolation, Mesh) {

    float float_max = std::numeric_limits<float>::max();

    std::vector<ThreeVector<float>> Vertex(4, ThreeVector<float>(0., 0., 0.));

    Vertex[0] = {0.0, 0.0, 0.0};
    Vertex[1] = {1.0, 0.0, 0.0};
    Vertex[2] = {0.0, 1.0, 0.0};
    Vertex[3] = {0.0, 0.0, 1.0};

    std::vector<ThreeVector<float>> DisplVector(4, ThreeVector<float>(0., 0., 0.));

    DisplVector[0] = {1.0, 2.0, 3.0};
    DisplVector[1] = {4.0, 5.0, 6.0};
    DisplVector[2] = {7.0, 8.0, 9.0};
    DisplVector[3] = {10.0, 11.0, 12.0};

    Laser LaserSet;

    LaserTrack Track1 = LaserTrack(Vertex,DisplVector);

    LaserSet.AppendTrack(Track1);

    ThreeVector<float> Location = {0., 0., 0.};

//    ASSERT_TRUE(Location == ThreeVector<float>(0.5, 0.5, 0.5));

// Create a array which contains the info of all 4 vertices of a cell
    std::array<std::pair<unsigned long, unsigned long>, 4> PointIndex;

// Initialize a displacement vector with zero
    ThreeVector<float> InterpolatedDispl = {0.0, 0.0, 0.0};

// Initialize Barycentric coordinate system (it will have 4 dimensions)
    std::vector<float> BaryCoord;

    Delaunay Mesh;
    Mesh = TrackMesher(LaserSet.GetTrackSet());

    // Find cell in the mesh where the point is located
    Delaunay::Cell_handle Cell = Mesh.locate(VectorToPoint(Location));

    // Loop over all four vertex points of the cell of interest
    for (unsigned vertex_no = 0; vertex_no < PointIndex.size(); vertex_no++) {
        // Get vertex info of the cell (track number, sample number)
        PointIndex[vertex_no] = Cell->vertex(vertex_no)->info();
    }

    // Initialize matrix for Location transformation into barycentric coordinate system
    Matrix3x3 TransMatrix = {{0, 0, 0},
                             {0, 0, 0},
                             {0, 0, 0}};

    // Loop over matrix rows
    for (unsigned row = 0; row < 3; row++) {
        // Loop over matrix columns
        for (unsigned column = 0; column < 3; column++) {
            // Fill transformation matrix elements
//            TransMatrix[row][column] = Vertex[column][row] - Vertex[3][row];
//            ASSERT_EQ(TransMatrix[row][column],)

            TransMatrix[row][column] =
                    LaserSet.GetTrackSet()[PointIndex[column].first].GetSamplePosition(PointIndex[column].second)[row] -
                    LaserSet.GetTrackSet()[PointIndex.back().first].GetSamplePosition(PointIndex.back().second)[row];

        }
    }

//    ASSERT_EQ(TransMatrix[0][0],0.);
//    ASSERT_EQ(TransMatrix[0][1],1.);
//    ASSERT_EQ(TransMatrix[0][2],0.);
//    ASSERT_EQ(TransMatrix[1][0],0.);
//    ASSERT_EQ(TransMatrix[1][1],0.);
//    ASSERT_EQ(TransMatrix[1][2],1.);
//    ASSERT_EQ(TransMatrix[2][0],-1.);
//    ASSERT_EQ(TransMatrix[2][1],-1.);
//    ASSERT_EQ(TransMatrix[2][2],-1.);




//    ThreeVector<float> RR4 = Location - Vertex[3];

    //// Reuse Location and store its position relative to the last vertex of the cell it is contained in
    //Location -= LaserMeshSet[PointIndex.back().first].GetSamplePosition(PointIndex.back().second);
    Location -= Vertex[3];

//    ASSERT_TRUE(Location == ThreeVector<float>(0.1, 0.1, -0.9));

    // If the transformation matrix can be successfully inverted
    if (TransMatrix.Invert()) {
        // Use inverted matrix to fill the first three coordinates
        ThreeVector<float> BC = TransMatrix * Location;
        BaryCoord = BC.GetStdVector();

//        ASSERT_NEAR(BC[0],1.,1E-3);
//        ASSERT_NEAR(BC[1],0.,1E-3);
//        ASSERT_NEAR(BC[2],0.,1E-3);



        // The sum of all barycentric coordinates has to be 1 by definition, use this to calculate the 4th coordinate
        BaryCoord.push_back(1 - BaryCoord[0] - BaryCoord[1] - BaryCoord[2]);

//        ASSERT_EQ(BC[0],0.1);
//        ASSERT_EQ(BC[1],-0.9);
//        ASSERT_EQ(BC[2],0.7);
//        ASSERT_EQ(BC[3],1.1);
    }
    else // if the matrix can't be inverted
    {
        //SAY SOMETHING IS WRONG!
        InterpolatedDispl = {-99.,-99.,-99.};

    }

//    ASSERT_NEAR(BaryCoord[0],1.,1E-3);
//    ASSERT_NEAR(BaryCoord[1],0.,1E-3);
//    ASSERT_NEAR(BaryCoord[2],0.,1E-3);
//    ASSERT_NEAR(BaryCoord[3],0.,1E-3);

    // Also barycentric coordinates need to be positive numbers (else the coordinate is outside of the cell).
    // So if one of the coordinates is smaller than zero
    if (BaryCoord[0] < 0.0 || BaryCoord[1] < 0.0 || BaryCoord[2] < 0.0 || BaryCoord[3] < 0.0) {
        InterpolatedDispl = {-999.,-999.,-999.};

    }

    // If the function is still alive, loop over all barycentric coordinates
    for (unsigned vertex_no = 0; vertex_no < 4; vertex_no++) {
        // Use the barycentric coordinates as a weight for the correction stored at this vertex in order to get the interpolated displacement
//        InterpolatedDispl += (LaserTrackSet[PointIndex[vertex_no].first].GetDisplacement(PointIndex[vertex_no].second) *
//        BaryCoord[vertex_no]);

        InterpolatedDispl += DisplVector[vertex_no] * BaryCoord[vertex_no];
    }

//    ASSERT_TRUE(InterpolatedDispl == ThreeVector<float>(1., 1., 1.));

//    ASSERT_NEAR(InterpolatedDispl[0],1.,1E-3);
//    ASSERT_NEAR(InterpolatedDispl[1],1.,1E-3);
//    ASSERT_NEAR(InterpolatedDispl[2],1.,1E-3);

//    ASSERT_EQ(InterpolatedDispl[0],1.);
    ASSERT_EQ(InterpolatedDispl[1],2.);
    ASSERT_EQ(InterpolatedDispl[2],3.);



}

int main(int ac, char* av[])
{
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}