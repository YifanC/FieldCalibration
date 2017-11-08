//
// Created by Yifan Chen on 09.10.17.
//

#include <limits>
#include "../include/ThreeVector.hpp"
#include "../include/LaserTrack.hpp"
#include "../include/Interpolation3D.hpp"
#include "../include/Laser.hpp"

std::pair<Laser, Laser>
DispLaserIteration(unsigned Nstep, Laser LaserSet1, Laser LaserSet2, bool CorrMapFlag){

    float float_max = std::numeric_limits<float>::max();
    Laser LaserRecoOrigin1 = LaserSet1;
    Laser LaserRecoOrigin2 = LaserSet2;

    for (int n = 0; n < Nstep; n++) {

        std::cout << "Processing correction step N " << n << " ... " << std::endl;

//        LaserSet1.CalcDisplacement(LaserTrack::ClosestPointCorr, Nstep - n);
//        LaserSet2.CalcDisplacement(LaserTrack::ClosestPointCorr, Nstep - n);

        // the direction of displacement is same as correction vector (reco to true)
        LaserSet1.CalcDisplacement(LaserTrack::ClosestPoint, Nstep - n);
        LaserSet2.CalcDisplacement(LaserTrack::ClosestPoint, Nstep - n);

        // At the last step, the biased track points should end on the true track lines
        if (n == (Nstep - 1)) {
            // the TRUE stands for opposite direction of the distortion direction (we calculate) and the correction direction (we will do here)
//            LaserSet1.AddCorrectionToReco(true);
//            LaserSet2.AddCorrectionToReco(true);

            // move the reco laser sets to the first step corrected laser tracks
            LaserSet1.AddCorrectionToReco();
            LaserSet2.AddCorrectionToReco();
        }
        else {

            Delaunay Mesh1 = TrackMesher(LaserSet1.GetTrackSet());
            Delaunay Mesh2 = TrackMesher(LaserSet2.GetTrackSet());

//            std::cout << "B" << std::difftime(std::time(NULL), timer) << " s" << std::endl;
//            std::cout << "total track " << LaserSet1.GetTrackSet().size() << std::endl;

            for (unsigned long track = 0; track < LaserSet1.GetTrackSet().size(); track++) {
//                std::cout << "Laser1:::Set--" << set << "--Nsetp--" << n << "--track--" << track << "--number--"
//                          << LaserSets1[set].GetTrackSet()[track].GetNumberOfSamples() << "||"
//                          << std::difftime(std::time(NULL), timer) << " s" << std::endl;
// reserve the space for the correction vector for each track
                unsigned long NrSamples1 = LaserSet1.GetTrackSet()[track].GetNumberOfSamples();
                std::vector<ThreeVector<float>> CorrPart1(NrSamples1, ThreeVector<float>(float_max, float_max, float_max));
//                        std::vector<ThreeVector<float>> CorrPart1(LaserSets1[set].GetTrackSet()[track].GetNumberOfSamples(),ThreeVector<float>(0,0,0));

// Loop over data points (samples) of each track
                // CorrPart is filled to 0 Threevector when the interpolation is failed !!!!! This is also a problem
                for (unsigned long sample = 0; sample < NrSamples1; sample++) {
                    CorrPart1[sample] = InterpolateCGAL(LaserSet2.GetTrackSet(), LaserSet2.GetTrackSet(), Mesh2,
                                                        LaserSet1.GetTrackSet()[track].GetSamplePosition(sample));
                }

                LaserSet1.GetTrackSet()[track].AddCorrectionToRecoPart(CorrPart1);
            }
//            std::cout << "F" << std::difftime(std::time(NULL), timer) << " s" << std::endl;

            for (unsigned long track = 0; track < LaserSet2.GetTrackSet().size(); track++) {
//                std::cout << "Laser2:::Set--" << set << "--Nsetp--" << Nstep << "--track--" << track
//                          << "--number--" << LaserSets2[set].GetTrackSet()[track].GetNumberOfSamples() << "||"
//                          << std::difftime(std::time(NULL), timer) << " s" << std::endl;
// reserve the space for the correction vector for each track
                unsigned long NrSamples2 = LaserSet2.GetTrackSet()[track].GetNumberOfSamples();
                std::vector<ThreeVector<float>> CorrPart2(NrSamples2, ThreeVector<float>(float_max, float_max, float_max));
//                        std::vector<ThreeVector<float>> CorrPart2(LaserSets2[set].GetTrackSet()[track].GetNumberOfSamples(),ThreeVector<float>(0,0,0));

// Loop over data points (samples) of each track
                for (unsigned long sample = 0; sample < NrSamples2; sample++) {
                    CorrPart2[sample] = InterpolateCGAL(LaserSet1.GetTrackSet(), LaserSet1.GetTrackSet(), Mesh1,
                                                        LaserSet2.GetTrackSet()[track].GetSamplePosition(sample));
                }
                LaserSet2.GetTrackSet()[track].AddCorrectionToRecoPart(CorrPart2);
            }
        }
    }
    LaserSet1.SetDisplacement(LaserRecoOrigin1, CorrMapFlag);
    LaserSet2.SetDisplacement(LaserRecoOrigin2, CorrMapFlag);

    std::pair<Laser, Laser> LaserWithDisplacement;
    LaserWithDisplacement = std::make_pair(LaserSet1,LaserSet2);

    return LaserWithDisplacement;

}

