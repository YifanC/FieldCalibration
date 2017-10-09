// C++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <chrono>
#include <cmath>
#include <cstring>
#include <thread>
#include <array>

// C headers
#include <pthread.h>
#include <unistd.h>
#include <getopt.h>

// ROOT headers
#include "TROOT.h"
#include "TObject.h"
#include "TApplication.h"
#include "TAttLine.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TProfile.h"
#include "TPolyLine3D.h"
#include "TStyle.h"
#include "TFrame.h"
#include "TFile.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TView3D.h"
#include "TTree.h"
#include "TChain.h"
#include "TVector3.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

// Own Files
#include "include/LaserTrack.hpp"
#include "include/ThreeVector.hpp"
#include "include/TPCVolumeHandler.hpp"
#include "include/Interpolation3D.hpp"
#include "include/Matrix3x3.hpp"
#include "include/Laser.hpp"
#include "include/Utilities.hpp"
#include "include/DriftVelocity.hpp"
#include "include/EfieldCalculation.hpp"

// Initialize functions defined below

//bool Twolasersys(int argc, char** argv);
Laser ReadRecoTracks(std::vector<std::string>);

void WriteRootFile(std::vector<ThreeVector<float>> &, TPCVolumeHandler &, std::string, bool CorrMapFlag);

void WriteTextFile(std::vector<ThreeVector<float>> &);

void LaserInterpThread(Laser &, const Laser &, const Delaunay &);

std::vector<Laser> ReachedExitPoint(const Laser &, float);

std::vector<ThreeVector<float>> Elocal(TPCVolumeHandler &, float cryoTemp, float E0, float v0, const char *);

std::vector<ThreeVector<float>> Eposition(TPCVolumeHandler &, float cryoTemp, float E0, float v0, const char *);

void WriteEmapRoot(std::vector<ThreeVector<float>> &Efield, TPCVolumeHandler &TPCVolume,
                   ThreeVector<unsigned long> Resolution, std::string);

// Set if the output displacement map is correction map (on reconstructed coordinate) or distortion map (on true coordinate)
// By default set it as correction map so we could continue calculate the E field map
bool CorrMapFlag = false; // Calculate Reco (coord) correction vectors for true; Calculate True (coord) distortion vectors for false
bool DoCorr = false; // Calculate Reco (coord) correction map for true; Skip calculation of True (coord) correction map for false
bool DoEmap = false; // Calculate electric map for true; Skip calculation of electric map for false
bool Merge2side = false;

// Main function
int main(int argc, char **argv) {
    // Start timer, just because it's nice to know how long shit takes
    time_t timer;
    std::time(&timer);

    // specify the amount of downsampling
    unsigned int n_split = 1;
    // Specify the number of steps for correction
    unsigned int Nstep = 1;

    // If there are to few input arguments, abort!
    if (argc < 2) {
        std::cerr << "ERROR: Too few arguments, use ./LaserCal <options> <input file names>" << std::endl;
        std::cerr << "options:  -d INTEGER  : Number of downsampling of the input dataset, default 1." << std::endl;
        return -1;
    }
    // Lets handle all options
    int c;
    while ((c = getopt(argc, argv, ":d:N:CDE")) != -1) {
        switch (c) {
            case 'd':
                n_split = atoi(optarg);
                break;
            case 'N':
                Nstep = atoi(optarg);
                break;
            case 'C':
                CorrMapFlag = true;
                break;
            case 'D':
                DoCorr = true;
                break;
            case 'E':
                DoEmap = true;
                break;
                // put in your case here. also add it to the while loop as an option or as required argument
        }
    }


//    // Now handle input files
//    std::vector<std::string> InputFiles;
//    unsigned int n_files = 0;
//    for (int i = optind; i < argc; i++) {
//        std::string filename(argv[i]);
//        // check if file exists
//        std::ifstream f(filename.c_str());
//        if (!f.good()) {
//            throw std::runtime_error(std::string("file does not exist: ") + filename);
//        }
//        InputFiles.push_back(filename);
//    }


    // Now handle input files
    std::vector<std::string> InputFiles1;
    std::vector<std::string> InputFiles2;
    unsigned int n_files = 0;
    for (int i = optind; i < argc; i++) {
        std::string filename(argv[i]);
        // check if file exists
        std::ifstream f(filename.c_str());
        if (!f.good()) {
            throw std::runtime_error(std::string("file does not exist: ") + filename);
        }

        TChain *tree = new TChain("lasers");
        tree->Add(filename.c_str());
        int side;
        tree->SetBranchAddress("side", &side);
        TCanvas *c1;
        tree->Draw("side>>hside", "");
        TH1F *hside = (TH1F *) gDirectory->Get("hside");
        int LCS = hside->GetMean();
//        c1->Close();
        delete tree;

        std::cout << "LCS: " << LCS << std::endl;

        if (LCS == 1) {
            InputFiles1.push_back(filename);
        } else if (LCS == 2) {
            InputFiles2.push_back(filename);
        } else {
            std::cerr << "The laser system is not labeled correctly." << std::endl;
        }
    }

    if (Merge2side) {
        InputFiles1.insert(InputFiles1.end(), InputFiles2.begin(), InputFiles2.end());
    } else {
        if (InputFiles1.empty() || InputFiles2.empty()) {
            std::cerr << "Please provide the laser data from 2 sides." << std::endl;
        }
    }

    // Choose detector dimensions, coordinate system offset and resolutions
    ThreeVector<float> DetectorSize = {256.04, 232.5, 1036.8};
    ThreeVector<float> DetectorOffset = {0.0, -DetectorSize[1] / static_cast<float>(2.0), 0.0};
    ThreeVector<unsigned long> DetectorResolution = {26, 26, 101};
    // Create the detector volume
    TPCVolumeHandler Detector(DetectorSize, DetectorOffset, DetectorResolution);

    ThreeVector<unsigned long> EMapResolution = {21, 21, 81};

    float cryoTemp = 89; // K
    float E0 = 0.273; // kV/cm
    float v0 = 1.11436; // mm/us, because of the fit of drift velocity as function of E field, while the LArSoft unit is cm/us

//    std::cout<<"The drift velocity range in consideration is from "<< ElectronDriftVelocity(cryoTemp, 0)<<" to "<< ElectronDriftVelocity(cryoTemp, 2*E0)<<std::endl;

    std::stringstream ss_outfile;
    std::stringstream ss_Eoutfile;
    float float_max = std::numeric_limits<float>::max();
    ThreeVector<float> Empty = {float_max, float_max, float_max};

//    // Set the name for Dmap
//    if (CorrMapFlag) {
//        ss_outfile << "RecoCorrection-N" << Nstep << "-S" << n_split << ".root";
//    }
//    if (!CorrMapFlag) {
//        ss_outfile << "TrueDistortion-N" << Nstep << "-S" << n_split << ".root";
//    }
//
//    ss_Eoutfile << "Emap-N" << Nstep << "-S" << n_split <<".root";

    ss_outfile << "RecoCorr-Simu.root";
    ss_Eoutfile << "Emap-Simu.root";

    if (DoCorr) {
        std::vector<std::vector<ThreeVector<float>>> DisplMapsHolder;
//        DisplMapsHolder.resize(n_split);

        float float_max = std::numeric_limits<float>::max();
        ThreeVector<float> Empty = {float_max, float_max, float_max};

        // Read data and store it to a Laser object
        std::cout << "Reading data..." << std::endl;
//        Laser FullTracks = ReadRecoTracks(InputFiles);
        Laser FullTracks1 = ReadRecoTracks(InputFiles1);
        Laser FullTracks2 = ReadRecoTracks(InputFiles2);

        // Here we split the laser set in multiple laser sets...
//        std::vector<Laser> LaserSets = SplitTrackSet(FullTracks, n_split);
        std::vector<Laser> LaserSets1 = SplitTrackSet(FullTracks1, n_split);
        std::vector<Laser> LaserSets2 = SplitTrackSet(FullTracks2, n_split);

        std::vector<Laser> LaserRecoOrigin1 = LaserSets1;
        std::vector<Laser> LaserRecoOrigin2 = LaserSets2;


        // Now we loop over each individual set and compute the displacement vectors.
        // TODO: This could be parallelized
        for (unsigned int set = 0; set < n_split; set++) {
            std::cout << "Processing subset " << set << " ... " << std::endl;

            // Calculate track displacement
            std::cout << "Find track displacements... " << std::endl;

            ////////////////////////////////////////////

//            Laser LaserRecoOrigin1 = LaserSets1;
//            Laser LaserRecoOrigin2 = LaserSets2;

            for (int n = 0; n < Nstep; n++) {

                std::cout << "Processing correction step N " << n << " ... " << std::endl;


                LaserSets1[set].CalcDisplacement(LaserTrack::ClosestPointCorr, Nstep - n);
                LaserSets2[set].CalcDisplacement(LaserTrack::ClosestPointCorr, Nstep - n);

                // when it becomes the last step, we require the biased track points will be dragged to the true track lines
                if (n == (Nstep - 1)) {
                    // the TRUE stands for opposite direction of the distortion direction (we calculate) and the correction direction (we will do here)
                    LaserSets1[set].AddCorrectionToReco(true);
                    LaserSets2[set].AddCorrectionToReco(true);
                } else {

                    std::cout << "A" << std::difftime(std::time(NULL), timer) << " s" << std::endl;

                    Delaunay Mesh1 = TrackMesher(LaserSets1[set].GetTrackSet());
                    Delaunay Mesh2 = TrackMesher(LaserSets2[set].GetTrackSet());

                    std::cout << "B" << std::difftime(std::time(NULL), timer) << " s" << std::endl;
                    std::cout << "total track " << LaserSets1[set].GetTrackSet().size() << std::endl;

                    for (unsigned long track = 0; track < LaserSets1[set].GetTrackSet().size(); track++) {
                        std::cout << "Laser1:::Set--" << set << "--Nsetp--" << n << "--track--" << track << "--number--"
                                  << LaserSets1[set].GetTrackSet()[track].GetNumberOfSamples() << "||"
                                  << std::difftime(std::time(NULL), timer) << " s" << std::endl;
                        // reserve the space for the correction vector for each track
                        std::vector<ThreeVector<float>> CorrPart1(
                                LaserSets1[set].GetTrackSet()[track].GetNumberOfSamples(),
                                ThreeVector<float>(float_max, float_max, float_max));
//                        std::vector<ThreeVector<float>> CorrPart1(LaserSets1[set].GetTrackSet()[track].GetNumberOfSamples(),ThreeVector<float>(0,0,0));

                        // Loop over data points (samples) of each track
                        for (unsigned long sample = 0;
                             sample < LaserSets1[set].GetTrackSet()[track].GetNumberOfSamples(); sample++) {
                            CorrPart1[sample] = InterpolateCGAL(LaserSets2[set].GetTrackSet(),
                                                                LaserSets2[set].GetTrackSet(), Mesh2,
                                                                LaserSets1[set].GetTrackSet()[track].GetSamplePosition(
                                                                        sample));
                        }

                        LaserSets1[set].GetTrackSet()[track].AddCorrectionToRecoPart(CorrPart1);
                    }
                    std::cout << "F" << std::difftime(std::time(NULL), timer) << " s" << std::endl;

                    for (unsigned long track = 0; track < LaserSets2[set].GetTrackSet().size(); track++) {
                        std::cout << "Laser2:::Set--" << set << "--Nsetp--" << Nstep << "--track--" << track
                                  << "--number--" << LaserSets2[set].GetTrackSet()[track].GetNumberOfSamples() << "||"
                                  << std::difftime(std::time(NULL), timer) << " s" << std::endl;
                        // reserve the space for the correction vector for each track
                        std::vector<ThreeVector<float>> CorrPart2(
                                LaserSets2[set].GetTrackSet()[track].GetNumberOfSamples(),
                                ThreeVector<float>(float_max, float_max, float_max));
//                        std::vector<ThreeVector<float>> CorrPart2(LaserSets2[set].GetTrackSet()[track].GetNumberOfSamples(),ThreeVector<float>(0,0,0));

                        // Loop over data points (samples) of each track
                        for (unsigned long sample = 0;
                             sample < LaserSets2[set].GetTrackSet()[track].GetNumberOfSamples(); sample++) {
                            CorrPart2[sample] = InterpolateCGAL(LaserSets1[set].GetTrackSet(),
                                                                LaserSets1[set].GetTrackSet(), Mesh1,
                                                                LaserSets2[set].GetTrackSet()[track].GetSamplePosition(
                                                                        sample));
                        }
                        LaserSets2[set].GetTrackSet()[track].AddCorrectionToRecoPart(CorrPart2);
                    }
                }
            }
            LaserSets1[set].SetDisplacement(LaserRecoOrigin1[set], CorrMapFlag);
            LaserSets2[set].SetDisplacement(LaserRecoOrigin2[set], CorrMapFlag);

            std::cout << "Time after N-step correction" << std::difftime(std::time(NULL), timer) << " s" << std::endl;

            ////////////////////////////////////////////

//            if (CorrMapFlag) {
//                // Suggestion: Choose ClosestPoint Algorithm
//                LaserSets[set].CalcDisplacement(LaserTrack::ClosestPointCorr);
//                // Now the laser data are based on the reconstructed coordinate.
//                // For CORRECTION MAP, no need to set the mesh on true space points
//            }
//
//            if (!CorrMapFlag) {
//                // Suggestion: Choose ClosestPoint Algorithm
//                LaserSets[set].CalcDisplacement(LaserTrack::ClosestPointDist);
//                // Now the laser tracks are based on the reconstructed coordinate.
//                // For DISTORTION MAP as output, set the mesh on the true space points
//                // the FALSE stands for opposite direction of the distortion direction (we calculate) and the correction direction (we will do here)
//                LaserSets[set].AddCorrectionToReco(false);
//            }

            // Create delaunay mesh
            std::cout << "Generate mesh..." << std::endl;
            Delaunay MeshMap1;
            Delaunay MeshMap2;

            std::cout << "Time after mesh " << std::difftime(std::time(NULL), timer) << " s" << std::endl;

            // The correction map is built on the mesh of reconstructed position which is the origin LaserSets
            if (CorrMapFlag) {
                MeshMap1 = TrackMesher(LaserRecoOrigin1[set].GetTrackSet());
                MeshMap2 = TrackMesher(LaserRecoOrigin2[set].GetTrackSet());
            }
                // The distortion map is built on the mesh of true position which is moved LaserSets
            else {
                MeshMap1 = TrackMesher(LaserSets1[set].GetTrackSet());
                MeshMap2 = TrackMesher(LaserSets2[set].GetTrackSet());
            }

//            Delaunay Mesh = TrackMesher(LaserSets[set].GetTrackSet());

            // Interpolate Displacement Map (regularly spaced grid)
            std::cout << "Start interpolation..." << std::endl;
            // LaserSets are now sitting on the true position, LaserRecoOrigin are sitting on the reco position

            // The correction map is based on reco space coord
            if (CorrMapFlag) {
                DisplMapsHolder.push_back(
                        InterpolateMap(LaserSets1[set].GetTrackSet(), LaserRecoOrigin1[set].GetTrackSet(), MeshMap1,
                                       Detector, CorrMapFlag));
                DisplMapsHolder.push_back(
                        InterpolateMap(LaserSets2[set].GetTrackSet(), LaserRecoOrigin2[set].GetTrackSet(), MeshMap2,
                                       Detector, CorrMapFlag));
            }
                // The distortion map is based on true space coord
            else {
                DisplMapsHolder.push_back(
                        InterpolateMap(LaserSets1[set].GetTrackSet(), LaserSets1[set].GetTrackSet(), MeshMap1, Detector,
                                       CorrMapFlag));
                DisplMapsHolder.push_back(
                        InterpolateMap(LaserSets2[set].GetTrackSet(), LaserSets2[set].GetTrackSet(), MeshMap2, Detector,
                                       CorrMapFlag));
            }
//            DisplMapsHolder[set] = InterpolateMap(LaserSets[set].GetTrackSet(), MeshMap1, Detector);
        }
        // Now we go on to create an unified displacement map
        std::vector<ThreeVector<float>> DisplacementMap(DisplMapsHolder.front().size(), ThreeVector<float>(0., 0., 0.));
        std::vector<float> Nvalid(DisplMapsHolder.front().size(), 0.);

        for (auto &SubMap: DisplMapsHolder) {
            for (unsigned int idx = 0; idx < DisplacementMap.size(); idx++) {
//                DisplacementMap[idx] = DisplacementMap[idx] + SubMap[idx] / static_cast<float>(n_split);
                if (SubMap[idx] != Empty) {
                    DisplacementMap[idx] = DisplacementMap[idx] + SubMap[idx];
                    Nvalid[idx]++;
                }
            }
        }

        for (unsigned int idx = 0; idx < DisplacementMap.size(); idx++) {
            if (Nvalid[idx] == 0) {
                // Set those bin with non valid number into float max again
                DisplacementMap[idx] = {float_max, float_max, float_max};
            } else {
                DisplacementMap[idx] = DisplacementMap[idx] / Nvalid[idx];
            }
        }

        // Fill displacement map into TH3 histograms and write them to file
        std::cout << "Write to File ..." << std::endl;
        WriteRootFile(DisplacementMap, Detector, ss_outfile.str(), CorrMapFlag);
    }

    // The Emap calculation works when the input is correction map
    if (DoEmap) {
        // The vector of Position and En must have the exactly the same index to make the interpolation (EInterpolateMap()) work
        auto E_field = Efield(Detector, cryoTemp, E0, v0, ss_outfile.str().c_str());
        std::vector<ThreeVector<float>> En = E_field.first;
        std::vector<ThreeVector<float>> Position = E_field.second;

//        std::vector<ThreeVector<float>> Position = Eposition(Detector, cryoTemp, E0, v0, ss_outfile.str().c_str());
//        std::vector<ThreeVector<float>> En = Elocal(Detector, cryoTemp, E0, v0, ss_outfile.str().c_str());

        // Create mesh for Emap
        std::cout << "Generate mesh for E field..." << std::endl;
        xDelaunay EMesh = Mesher(Position, Detector);

        // Interpolate E Map (regularly spaced grid)
        std::cout << "Start interpolation the E field..." << std::endl;
        std::vector<ThreeVector<float>> EMap = EInterpolateMap(En, Position, EMesh, Detector, EMapResolution);

        // Fill displacement map into TH3 histograms and write them to file
        std::cout << "Write Emap to File ..." << std::endl;
        WriteEmapRoot(EMap, Detector, EMapResolution, ss_Eoutfile.str());
    }


    std::cout << "End of program after " << std::difftime(std::time(NULL), timer) << " s" << std::endl;


} // end main


Laser ReadRecoTracks(std::vector<std::string> InputFiles) {
    // Create Laser (collection of laser tracks) this will be the returned object
    Laser TrackSelection;

    // Initialize read variables, the pointers for more complex data structures 
    // are very important for Root. Rene Brun in hell (do you see what I did there?)
    int EventNumber;
    std::vector<TVector3> TrackSamples;
    std::vector<TVector3> *pTrackSamples = &TrackSamples;
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


void WriteRootFile(std::vector<ThreeVector<float>> &InterpolationData, TPCVolumeHandler &TPCVolume,
                   std::string OutputFilename, bool CorrMapFlag = false) {
    // Store TPC properties which are important for the TH3 generation

    ThreeVector<unsigned long> Resolution = TPCVolume.GetDetectorResolution();
    ThreeVector<float> MinimumCoord = TPCVolume.GetMapMinimum();
    ThreeVector<float> MaximumCoord = TPCVolume.GetMapMaximum();
    ThreeVector<float> Unit = {TPCVolume.GetDetectorSize()[0] / (Resolution[0] - 1),
                               TPCVolume.GetDetectorSize()[1] / (Resolution[1] - 1),
                               TPCVolume.GetDetectorSize()[2] / (Resolution[2] - 1)};


//    std::cout<<"MinX: "<<MinimumCoord[0]<<"; MinY: "<<MinimumCoord[1]<<"; MinZ: "<<MinimumCoord[2]<<std::endl;
//    std::cout<<"MaxX: "<<MaximumCoord[0]<<"; MaxY: "<<MaximumCoord[1]<<"; MaxZ: "<<MaximumCoord[2]<<std::endl;

    // For Reco coord based
    unsigned Extension = 0;
    if (CorrMapFlag) {
        Extension = 2; //Must be consistent with the setup as "InterpolationData"
    }
    // Initialize all TH3F
    std::vector<TH3F> RecoDisplacement;
    //2 in "2*Extension" stands for extension of both sides
    RecoDisplacement.push_back(TH3F("Reco_Displacement_X", "Reco Displacement X", Resolution[0] + 2 * Extension,
                                    MinimumCoord[0] - Unit[0] * (Extension + 0.5),
                                    MaximumCoord[0] + Unit[0] * (Extension + 0.5), Resolution[1] + 2 * Extension,
                                    MinimumCoord[1] - Unit[1] * (Extension + 0.5),
                                    MaximumCoord[1] + Unit[1] * (Extension + 0.5), Resolution[2] + 2 * Extension,
                                    MinimumCoord[2] - Unit[2] * (Extension + 0.5),
                                    MaximumCoord[2] + Unit[2] * (Extension + 0.5)));
    RecoDisplacement.push_back(TH3F("Reco_Displacement_Y", "Reco Displacement Y", Resolution[0] + 2 * Extension,
                                    MinimumCoord[0] - Unit[0] * (Extension + 0.5),
                                    MaximumCoord[0] + Unit[0] * (Extension + 0.5), Resolution[1] + 2 * Extension,
                                    MinimumCoord[1] - Unit[1] * (Extension + 0.5),
                                    MaximumCoord[1] + Unit[1] * (Extension + 0.5), Resolution[2] + 2 * Extension,
                                    MinimumCoord[2] - Unit[2] * (Extension + 0.5),
                                    MaximumCoord[2] + Unit[2] * (Extension + 0.5)));
    RecoDisplacement.push_back(TH3F("Reco_Displacement_Z", "Reco Displacement Z", Resolution[0] + 2 * Extension,
                                    MinimumCoord[0] - Unit[0] * (Extension + 0.5),
                                    MaximumCoord[0] + Unit[0] * (Extension + 0.5), Resolution[1] + 2 * Extension,
                                    MinimumCoord[1] - Unit[1] * (Extension + 0.5),
                                    MaximumCoord[1] + Unit[1] * (Extension + 0.5), Resolution[2] + 2 * Extension,
                                    MinimumCoord[2] - Unit[2] * (Extension + 0.5),
                                    MaximumCoord[2] + Unit[2] * (Extension + 0.5)));



    // Loop over all xbins
    for (unsigned xbin = 0; xbin < Resolution[0] + 2 * Extension; xbin++) {
        // Loop over all ybins
        for (unsigned ybin = 0; ybin < Resolution[1] + 2 * Extension; ybin++) {
            // Loop over all zbins
            for (unsigned zbin = 0; zbin < Resolution[2] + 2 * Extension; zbin++) {
                // Loop over all coordinates

                for (unsigned coord = 0; coord < 3; coord++) {
                    // Fill interpolated grid points into histograms
                    // bin=0 is underflow, bin = nbin+1 is overflow
                    RecoDisplacement[coord].SetBinContent(xbin + 1, ybin + 1, zbin + 1, InterpolationData[zbin +
                                                                                                          (Resolution[2] *
                                                                                                           (ybin +
                                                                                                            Resolution[1] *
                                                                                                            xbin))][coord]);
                    // It's equivalent to the following expression
                    // Remember, the range of the hist bin is (1, nbins), while when we fill the vector, it starts from 0. (0,nbins-1)
                    // RecoDisplacement[coord].SetBinContent(xbin+1,ybin+1,zbin+1, InterpolationData[zbin+ybin*TPCVolume.GetDetectorResolution()[2]+xbin*TPCVolume.GetDetectorResolution()[1]*TPCVolume.GetDetectorResolution()[2]][coord]);
                } // end coordinate loop
            } // end zbin loop
        } // end ybin loop
    } // end zbin loop

    // Open and recreate output file

    TFile OutputFile(OutputFilename.c_str(), "recreate");

    // Loop over space coordinates
    for (unsigned coord = 0; coord < RecoDisplacement.size(); coord++) {
        // Write every TH3 map into file
        RecoDisplacement[coord].Write();
    }

    // Close output file and clean up
    OutputFile.Close();
    gDirectory->GetList()->Delete();
}


void WriteTextFile(std::vector<ThreeVector<float>> &InterpolationData) {
    // Initialize stream to file
    std::ofstream OutputFile;

    // Open output file
    OutputFile.open("Reco.txt", std::ios::out);

    // Loop over all interpolated data points
    for (unsigned entry = 0; entry < InterpolationData.size(); entry++) {
        // Write every point into a seperate line
        OutputFile << InterpolationData[entry][0] << InterpolationData[entry][1] << InterpolationData[entry][2];
    }

    // Close file
    OutputFile.close();
} // WriteRootFile


// This is the multi-threading interpolation function. Just hangs out here, for legacy purposes 
void LaserInterpThread(Laser &LaserTrackSet, const Laser &InterpolationLaser, const Delaunay &InterpolationMesh) {
    LaserTrackSet.InterpolateTrackSet(InterpolationLaser, InterpolationMesh);
} // LaserInterpThread

// Split the laser track set into tracks that reached the expected exit point (within a configurable region) and others.
// First entry of the return vector is tracks that reach the exit point, second is the ones that do not reach it.


//// The root file does not have to be the argument
//std::vector<ThreeVector<float>>
//Elocal(TPCVolumeHandler &TPCVolume, float cryoTemp, float E0, float v0, const char *root_name) {
////    TFile *InFile = new TFile("RecoCorrection.root","READ");
//    TFile *InFile = new TFile(root_name, "READ");
//
//    TH3F *Dx = (TH3F *) InFile->Get("Reco_Displacement_X");
//    TH3F *Dy = (TH3F *) InFile->Get("Reco_Displacement_Y");
//    TH3F *Dz = (TH3F *) InFile->Get("Reco_Displacement_Z");
//
//    ThreeVector<unsigned long> Resolution = TPCVolume.GetDetectorResolution();
//    ThreeVector<float> DetectorReso = {TPCVolume.GetDetectorSize()[0] / (Resolution[0] - 1),
//                                       TPCVolume.GetDetectorSize()[1] / (Resolution[1] - 1),
//                                       TPCVolume.GetDetectorSize()[2] / (Resolution[2] - 1)};
//    float Delta_x = DetectorReso[0]; //cm
//
////    std::vector<ThreeVector< float>> En(TPCVolume.GetDetectorResolution()[2] * TPCVolume.GetDetectorResolution()[1] * (TPCVolume.GetDetectorResolution()[0]-1));
//    std::vector<ThreeVector<float>> En;
//
//    for (unsigned zbin = 0; zbin < TPCVolume.GetDetectorResolution()[2]; zbin++) {
//        for (unsigned ybin = 0; ybin < TPCVolume.GetDetectorResolution()[1]; ybin++) {
//            // the number of x bin in Emap is one less than in the displacement map
//            // because we only consider the gap here
//            for (unsigned xbin = 0; xbin < (TPCVolume.GetDetectorResolution()[0] - 1); xbin++) {
//                //xbin =1, x=5 close to the anode; xbin = Nx, x=250 close to the cathode
//                // the corner has the weight of the bin. Do not move the weight to the geometry center!
//                ThreeVector<float> RecoGrid(xbin * DetectorReso[0],
//                                            ybin * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
//                                            zbin * DetectorReso[2]);
//                ThreeVector<float> Dxyz = {(float) Dx->GetBinContent(xbin + 1, ybin + 1, zbin + 1),
//                                           (float) Dy->GetBinContent(xbin + 1, ybin + 1, zbin + 1),
//                                           (float) Dz->GetBinContent(xbin + 1, ybin + 1, zbin + 1)};
//                ThreeVector<float> True = RecoGrid + Dxyz;
//
//                ThreeVector<float> RecoGrid_next((xbin + 1) * DetectorReso[0],
//                                                 ybin * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
//                                                 zbin * DetectorReso[2]);
//                ThreeVector<float> Dxyz_next = {(float) Dx->GetBinContent(xbin + 2, ybin + 1, zbin + 1),
//                                                (float) Dy->GetBinContent(xbin + 2, ybin + 1, zbin + 1),
//                                                (float) Dz->GetBinContent(xbin + 2, ybin + 1, zbin + 1)};
//                ThreeVector<float> True_next = RecoGrid_next + Dxyz_next;
//
//                ThreeVector<float> Rn = True_next - True;
//                float vn = Rn.GetNorm() / Delta_x * v0; // mm/us, the magnitude of the drift velocity at the local point
//                //To be very careful that Rn.GetNorm() and Delta_x are in the unit of cm, not mm
//
////                if(searchE(vn,cryoTemp,E0)>0.6){
////                    std::cout<<"Need investigation! xbin: "<<xbin<<"; ybin: "<<ybin<<"; zbin: "<<zbin<<"; |Rn|: "<<Rn.GetNorm()<<"; |E|: "<<searchE(vn,cryoTemp,E0)<<std::endl;
////                }
//
//                // the E field as a vector has the same direction of Rn (Threevector)
//                if (searchE(vn, cryoTemp, E0) < 0.5 * std::numeric_limits<float>::max()) {
//                    En.push_back(searchE(vn, cryoTemp, E0) / Rn.GetNorm() * Rn);
//                }
//
////                En[xbin+ybin*(Resolution[0]-1)+zbin*(Resolution[0]-1)*Resolution[1]] = searchE(vn,cryoTemp,E0) / Rn.GetNorm() * Rn;
//            }
//        }
//    }
//
//    std::cout << "En size: " << En.size() << std::endl;
//    return En;
//}
//
//std::vector<ThreeVector<float>>
//Eposition(TPCVolumeHandler &TPCVolume, float cryoTemp, float E0, float v0, const char *root_name) {
//    ThreeVector<unsigned long> Resolution = TPCVolume.GetDetectorResolution();
//    ThreeVector<float> DetectorReso = {TPCVolume.GetDetectorSize()[0] / (Resolution[0] - 1),
//                                       TPCVolume.GetDetectorSize()[1] / (Resolution[1] - 1),
//                                       TPCVolume.GetDetectorSize()[2] / (Resolution[2] - 1)};
//
//    std::cout << "name: " << root_name << std::endl;
//    TFile *InFile = new TFile(root_name, "READ");
//
//    TH3F *Dx = (TH3F *) InFile->Get("Reco_Displacement_X");
//    TH3F *Dy = (TH3F *) InFile->Get("Reco_Displacement_Y");
//    TH3F *Dz = (TH3F *) InFile->Get("Reco_Displacement_Z");
//    float Delta_x = DetectorReso[0]; //cm
//
////    std::vector<ThreeVector<float>> Position(Resolution[2]*Resolution[1]*(Resolution[0]-1));
//    std::vector<ThreeVector<float>> Position;
//
//    // the position should be consistent to the one in the EInterpolateMap()
//    for (unsigned zbin = 0; zbin < Resolution[2]; zbin++) {
//        for (unsigned ybin = 0; ybin < Resolution[1]; ybin++) {
//            // since we calculate Elocal by the gap, the number of x bin in Emap is one less than in the displacement map
//            for (unsigned xbin = 0; xbin < (Resolution[0] - 1); xbin++) {
//                //xbin =1, x=5 close to the anode; xbin = Nx, x=250 close to the cathode
//                ThreeVector<float> RecoGrid(xbin * DetectorReso[0],
//                                            ybin * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
//                                            zbin * DetectorReso[2]);
//                ThreeVector<float> Dxyz = {(float) Dx->GetBinContent(xbin + 1, ybin + 1, zbin + 1),
//                                           (float) Dy->GetBinContent(xbin + 1, ybin + 1, zbin + 1),
//                                           (float) Dz->GetBinContent(xbin + 1, ybin + 1, zbin + 1)};
//                ThreeVector<float> True = RecoGrid + Dxyz;
//
//                ThreeVector<float> RecoGrid_next((xbin + 1) * DetectorReso[0],
//                                                 ybin * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
//                                                 zbin * DetectorReso[2]);
//                ThreeVector<float> Dxyz_next = {(float) Dx->GetBinContent(xbin + 2, ybin + 1, zbin + 1),
//                                                (float) Dy->GetBinContent(xbin + 2, ybin + 1, zbin + 1),
//                                                (float) Dz->GetBinContent(xbin + 2, ybin + 1, zbin + 1)};
//                ThreeVector<float> True_next = RecoGrid_next + Dxyz_next;
//
//                ThreeVector<float> Rn = True_next - True;
//
//                // To synchronize the vector order of En and Position through the output of searchE (if the output is floatmax, then abandon both elements in Position and En)
//                float vn = Rn.GetNorm() / Delta_x * v0;
//                if (searchE(vn, cryoTemp, E0) < 0.5 * std::numeric_limits<float>::max()) {
//                    Position.push_back(searchE(vn, cryoTemp, E0) / Rn.GetNorm() * Rn);
//                }
////                Position[xbin+ybin*(Resolution[0]-1)+zbin*(Resolution[0]-1)*Resolution[1]] = True + (float) 0.5 * Rn;
//            }
//        }
//    }
//    std::cout << "Position size: " << Position.size() << std::endl;
//    return Position;
//}

// Write Emap into TH3 and store in root file
void WriteEmapRoot(std::vector<ThreeVector<float>> &Efield, TPCVolumeHandler &TPCVolume,
                   ThreeVector<unsigned long> Resolution, std::string OutputFilename) {
    // Store TPC properties which are important for the TH3 generation
//    ThreeVector<unsigned long> Resolution = {21,21,81};
//    ThreeVector<unsigned long> Resolution = TPCVolume.GetDetectorResolution();
    ThreeVector<float> MinimumCoord = TPCVolume.GetMapMinimum();
    ThreeVector<float> MaximumCoord = TPCVolume.GetMapMaximum();
    ThreeVector<float> Unit = {TPCVolume.GetDetectorSize()[0] / (Resolution[0] - 1),
                               TPCVolume.GetDetectorSize()[1] / (Resolution[1] - 1),
                               TPCVolume.GetDetectorSize()[2] / (Resolution[2] - 1)};

    // Initialize all TH3F
    std::vector<TH3F> Emap;
    Emap.push_back(TH3F("Emap_X", "E field map X", Resolution[0], MinimumCoord[0] - Unit[0] * 0.5,
                        MaximumCoord[0] + Unit[0] * 0.5, Resolution[1], MinimumCoord[1] - Unit[1] * 0.5,
                        MaximumCoord[1] + Unit[1] * 0.5, Resolution[2], MinimumCoord[2] - Unit[2] * 0.5,
                        MaximumCoord[2] + Unit[2] * 0.5));
    Emap.push_back(TH3F("Emap_Y", "E field map Y", Resolution[0], MinimumCoord[0] - Unit[0] * 0.5,
                        MaximumCoord[0] + Unit[0] * 0.5, Resolution[1], MinimumCoord[1] - Unit[1] * 0.5,
                        MaximumCoord[1] + Unit[1] * 0.5, Resolution[2], MinimumCoord[2] - Unit[2] * 0.5,
                        MaximumCoord[2] + Unit[2] * 0.5));
    Emap.push_back(TH3F("Emap_Z", "E field map Z", Resolution[0], MinimumCoord[0] - Unit[0] * 0.5,
                        MaximumCoord[0] + Unit[0] * 0.5, Resolution[1], MinimumCoord[1] - Unit[1] * 0.5,
                        MaximumCoord[1] + Unit[1] * 0.5, Resolution[2], MinimumCoord[2] - Unit[2] * 0.5,
                        MaximumCoord[2] + Unit[2] * 0.5));

    // the loop should be consistent to the one in the EInterpolateMap()
    for (unsigned xbin = 0; xbin < Resolution[0]; xbin++) {
        for (unsigned ybin = 0; ybin < Resolution[1]; ybin++) {
            for (unsigned zbin = 0; zbin < Resolution[2]; zbin++) {
                // Loop over all coordinates dx,dy,dz
                for (unsigned coord = 0; coord < 3; coord++) {
                    // Fill interpolated grid points into histograms. bin=0 is underflow, bin = nbin+1 is overflow
                    Emap[coord].SetBinContent(xbin + 1, ybin + 1, zbin + 1, Efield[zbin + ybin * Resolution[2] +
                                                                                   xbin * Resolution[2] *
                                                                                   Resolution[1]][coord]);
                } // end coordinate loop
//                std::cout<<"xbin: "<<xbin<<"; ybin: "<<ybin<<"; zbin: "<<zbin<<"---Ex: "<<Efield[zbin+ybin*Resolution[2]+xbin*Resolution[2]*Resolution[1]][0]<<"; Ey: "<<Efield[zbin+ybin*Resolution[2]+xbin*Resolution[2]*Resolution[1]][1]<<"; Ez: "<< Efield[zbin+ybin*Resolution[2]+xbin*Resolution[2]*Resolution[1]][2]<<std::endl;
            } // end zbin loop
        } // end ybin loop
    } // end zbin loop

    // Open and recreate output file
    TFile OutputFile(OutputFilename.c_str(), "recreate");

    // Loop over space coordinates
    for (unsigned coord = 0; coord < Emap.size(); coord++) {
        // Write every TH3 map into file
        Emap[coord].Write();
    }

    // Close output file and clean up
    OutputFile.Close();
    gDirectory->GetList()->Delete();
}