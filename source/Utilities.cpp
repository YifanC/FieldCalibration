//
// Created by matthias on 28.04.17.
//

#include "../include/Utilities.hpp"
#include "../include/Laser.hpp"
#include "TROOT.h"
#include "TChain.h"
#include <fstream>
#include <iostream>

std::vector<Laser> ReachedExitPoint(const Laser &LaserSet, float ExitBoundary) {
    /*
     * Function that splits the input data set into two sets. The first set contains the
     * tracks that have a distance between the last reconstructed point and the expected
     * exit point smaller than the specified distance (ExitBoundary). All other LaserTracks
     * will be put in the second set.
     */

    std::vector<Laser> Selection;
    Selection.resize(2);

    for (auto &Track : LaserSet.GetTrackSet()) {
        auto ExitToLast = Track.GetExitPoint() - Track.GetBack();
        auto d = ExitToLast.GetNorm();

        if (d < ExitBoundary) {
            Selection.front().AppendTrack(Track);
        } else {
            Selection.back().AppendTrack(Track);
        }
    }
    return Selection;
}

std::vector<Laser> SplitTrackSet(const Laser &LaserSet, unsigned int Downsample) {
    /*
     * This function separates the input laser data set into multiple smaller data sets,
     * the number of produced set is specified by the Downsample value.
     */
    std::vector<Laser> Sets;
    Sets.resize(Downsample);

    for (auto &Track : LaserSet.GetTrackSet()) {

        auto SourceTrack = Track.GetReco();

        for (unsigned long offset = 0; offset < Downsample; offset++) {
            std::vector<ThreeVector<float>> SampledRecoTrack;

            for (unsigned long idx = offset; idx < SourceTrack.size(); idx += Downsample) {
                SampledRecoTrack.push_back(SourceTrack[idx]);
            }
            LaserTrack SampledTrack(Track.GetEntryPoint(), Track.GetExitPoint(), SampledRecoTrack);
            Sets[offset].AppendTrack(SampledTrack);
            SampledRecoTrack.clear();

        }
    }
    return Sets;
}

std::vector<Laser> InterlacedIteration(const Laser &LaserSet) {
    /* This function separates the input laser data set into 2 data subsamples for interlaced iterated correction.*/
    std::vector<Laser> Sets;
    Sets.resize(2);

    unsigned int number = 0;

    for (auto &Track : LaserSet.GetTrackSet()) {

        if(number%2 == 0){
            Sets[0].AppendTrack(Track);
        }
        if(number%2 == 1){
            Sets[1].AppendTrack(Track);
        }

        number++;
    }
    return Sets;
}

Laser ReadRecoTracks(std::vector<std::string> InputFiles) {
    // Create Laser (collection of laser tracks) this will be the returned object
    Laser TrackSelection;

    // Initialize read variables, the pointers for more complex data structures
    // are very important for Root. Rene Brun in hell (do you see what I did there?)
    int EventNumber;

    std::vector<TVector3> TrackSamples;
    std::vector<TVector3>* pTrackSamples = &TrackSamples;

    TVector3 EntryPoint;
    TVector3 *pEntryPoint = &EntryPoint;
    TVector3 ExitPoint;
    TVector3 *pExitPoint = &ExitPoint;

    // Open TChains to store all trees
    TChain *LaserInfoTree = new TChain("lasers");
    TChain *RecoTrackTree = new TChain("tracks");

    // Loop through all input files and add them to the TChain
    for (auto const &InFile : InputFiles) {
        // Open input file and add to TChains
        LaserInfoTree->Add(InFile.c_str());
        RecoTrackTree->Add(InFile.c_str());
    }

    // Assign branch addresses
    LaserInfoTree->SetBranchAddress("entry", &pEntryPoint);
    LaserInfoTree->SetBranchAddress("exit", &pExitPoint);
    RecoTrackTree->SetBranchAddress("track", &pTrackSamples);
    RecoTrackTree->SetBranchAddress("event", &EventNumber);

    // Only start read out when both trees have the same amount of entries
    if (LaserInfoTree->GetEntries() == RecoTrackTree->GetEntries()) {
        // Loop over all tree entries
        for (Size_t tree_index = 0; tree_index < RecoTrackTree->GetEntries(); tree_index++) {
            // Get tree entries of both trees
            LaserInfoTree->GetEntry(tree_index);
            RecoTrackTree->GetEntry(tree_index);

            // Sorting wouldn't change the track physically
            // For closestpoint method, it doesn't matter, while to derivative method yes
            // But for the moment, when reconstruction has a big problem, it is not encouraged to use derivative method

            // This here sorts the tracks by their distance to the EntryPoint. The algorithm uses a lambda
            // It will compare the distance to the EntryPoint of two vector entries A & B
            std::sort(TrackSamples.begin(), TrackSamples.end(), [&EntryPoint](TVector3 A, TVector3 B) {
                          A -= EntryPoint;
                          B -= EntryPoint;
                          // Here only the squared distance was used to avoid costly sqrt operations
                          return A.Mag2() > B.Mag2();
                      }
            );

            // This step will erase all double entries. First std::unique shifts every double to the end
            // of the vector and gives back the new end point of the data set. After that we erase the HistRange
            // between this new end and the real end of the vector
            TrackSamples.erase(std::unique(TrackSamples.begin(), TrackSamples.end()), TrackSamples.end());

            // Add new track to Laser TrackSelection
            TrackSelection.AppendTrack(LaserTrack(EntryPoint, ExitPoint, TrackSamples));
        }
    } else // If the trees don't have the same amount of entries, through error (I know not propper error handling)
    {
        std::cerr << "ERROR: Two TTrees don't have the same amount of entries!" << std::endl;
    }

//     delete pEntryPoint;
//     delete pExitPoint;
//     delete pTrackSamples;

    gDirectory->GetList()->Delete();

    delete LaserInfoTree;
    delete RecoTrackTree;

    return TrackSelection;
} // end ReadRecoTracks


// Merge two "Laser"s into one "Laser"
// the merged Laser will be independent to LaserA and LaserB
Laser MergeLaser(const Laser &LaserA, const Laser &LaserB) {

    // Initialize output of type Laser
    Laser MergedLaser;

    MergedLaser.Mergewith(LaserA);
    MergedLaser.Mergewith(LaserB);

    return MergedLaser;
}