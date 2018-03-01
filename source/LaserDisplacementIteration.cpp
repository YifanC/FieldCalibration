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

        // the direction of displacement is same as correction vector (reco to true)
        LaserSet1.CalcDisplacement(LaserTrack::ClosestPoint, Nstep - n);
        LaserSet2.CalcDisplacement(LaserTrack::ClosestPoint, Nstep - n);

        // At the last step, the biased track points should end on the true track lines
        if (n == (Nstep - 1)) {

            // move the reco laser sets to the first step corrected laser tracks
            LaserSet1.AddCorrectionToReco();
            LaserSet2.AddCorrectionToReco();
        }
        else {

            Delaunay Mesh1 = TrackMesher(LaserSet1.GetTrackSet());
            Delaunay Mesh2 = TrackMesher(LaserSet2.GetTrackSet());

            // LaserSet1 is going to change during iteration, while we need "origin" LaserSet1 of this step for Laser2 Mesh
            // So the order of LaserSet1 or 2 matters here
            auto LaserSet1_stepCopy = LaserSet1.GetTrackSet();

            for (unsigned long track = 0; track < LaserSet1.GetTrackSet().size(); track++) {

                // reserve the space for the correction vector for each track
                unsigned long NrSamples1 = LaserSet1.GetTrackSet()[track].GetNumberOfSamples();
                std::vector<ThreeVector<float>> CorrPart1(NrSamples1, ThreeVector<float>(float_max, float_max, float_max));

                // Loop over data points (samples) of each track
                // TODO CorrPart is filled to 0 Threevector when the interpolation is failed !!!!! This is also a problem

                for (unsigned long sample = 0; sample < NrSamples1; sample++) {
//                    CorrPart1[sample] = InterpolateCGAL(LaserSet2.GetTrackSet(), LaserSet2.GetTrackSet(), Mesh2,
//                                                        LaserSet1.GetTrackSet()[track].GetSamplePosition(sample));
                    CorrPart1[sample] = InterpolateCGAL(LaserSet2.GetTrackSet(), Mesh2,
                                                        LaserSet1.GetTrackSet()[track].GetSamplePosition(sample));
                }

                LaserSet1.GetTrackSet()[track].AddCorrectionToRecoPart(CorrPart1);
            }

            for (unsigned long track = 0; track < LaserSet2.GetTrackSet().size(); track++) {

                // reserve the space for the correction vector for each track
                unsigned long NrSamples2 = LaserSet2.GetTrackSet()[track].GetNumberOfSamples();
                std::vector<ThreeVector<float>> CorrPart2(NrSamples2, ThreeVector<float>(float_max, float_max, float_max));

                // Loop over data points (samples) of each track
                for (unsigned long sample = 0; sample < NrSamples2; sample++) {
//                    CorrPart2[sample] = InterpolateCGAL(LaserSet1.GetTrackSet(), LaserSet1.GetTrackSet(), Mesh1,
//                                                        LaserSet2.GetTrackSet()[track].GetSamplePosition(sample));
                    CorrPart2[sample] = InterpolateCGAL(LaserSet1_stepCopy, Mesh1,
                                                        LaserSet2.GetTrackSet()[track].GetSamplePosition(sample));
                }
                LaserSet2.GetTrackSet()[track].AddCorrectionToRecoPart(CorrPart2);
            }
        }
    }

    std::pair<Laser, Laser> LaserWithDisplacement;

    if(CorrMapFlag){
        // Set correction vectors on reconstructed track for distortion map
        LaserRecoOrigin1.SetDisplacement(LaserSet1);
        LaserRecoOrigin2.SetDisplacement(LaserSet2);
        LaserWithDisplacement = std::make_pair(LaserRecoOrigin1,LaserRecoOrigin2);

    } else{
        // Set distortion vectors on true track for distortion map
        LaserSet1.SetDisplacement(LaserRecoOrigin1);
        LaserSet2.SetDisplacement(LaserRecoOrigin2);
        LaserWithDisplacement = std::make_pair(LaserSet1,LaserSet2);
    }

    return LaserWithDisplacement;


//    std::pair<Laser, Laser> LaserWithDisplacement;
//    LaserWithDisplacement = std::make_pair(LaserSet1,LaserSet2);
//
//    return LaserWithDisplacement;

//    LaserSet1.SetDisplacement(LaserRecoOrigin1, CorrMapFlag);
//    LaserSet2.SetDisplacement(LaserRecoOrigin2, CorrMapFlag);
//
//    std::pair<Laser, Laser> LaserWithDisplacement;
//    LaserWithDisplacement = std::make_pair(LaserSet1,LaserSet2);
//
//    return LaserWithDisplacement;

}

