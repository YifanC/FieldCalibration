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

/*
#ifdef _OPENMP
#include "omp.h"
#endif
*/

#include <dirent.h>

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
#include "include/LaserDisplacementIteration.hpp"
#include "include/WeightAverage.hpp"

// Initialize functions defined below

Laser ReadRecoTracks(std::vector<std::string>);

LaserTrack Anode(TPCVolumeHandler &TPCVolume);

void WriteRootFile(std::vector<ThreeVector<float>> &, TPCVolumeHandler &, std::string);

void WriteRootFileDeviation(std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>> &, TPCVolumeHandler &, std::string);

void WriteTextFile(std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>> &, std::string);

void LaserInterpThread(Laser &, const Laser &, const Delaunay &);

std::vector<Laser> ReachedExitPoint(const Laser &, float);

void WriteEmapRoot(std::vector<ThreeVector<float>> &Efield, TPCVolumeHandler &TPCVolume,
                   ThreeVector<unsigned long> Resolution, float E0, std::string);

// Set if the output displacement map is correction map (on reconstructed coordinate) or distortion map (on true coordinate)
// By default set it as correction map so we could continue calculate the E field map
bool CorrMapFlag = false; // Calculate Reco (coord) correction vectors for true; Calculate True (coord) distortion vectors for false
bool DoCorr = false; // Calculate Reco (coord) correction map for true; Skip calculation of True (coord) correction map for false
bool DoEmap = false; // Calculate electric map for true; Skip calculation of electric map for false
bool TwoSideIter = false;
bool InterlacedIter = false;
bool CosmicLaserIter = false;
bool DBoundary = false;
bool EBoundary = false;
bool WeightAverage = false;

// Main function
int main(int argc, char **argv) {

    // Start timer, just because it's nice to know how long shit takes
    time_t timer;
    std::time(&timer);

    // specify the amount of downsampling
    unsigned int n_split = 1;
    unsigned int n_threads = 1;
    // Specify the number of iteration steps. If Nstep = 1, there will be no iteration.
    unsigned int Nstep = 1;

    // If there are to few input arguments, abort!
    if (argc < 2) {
        std::cerr << "ERROR: Too few arguments, use ./LaserCal <options> <input file names>" << std::endl;
        std::cerr << "options:  -d INTEGER  : Number of downsampling of the input dataset, default 1." << std::endl;
        std::cerr << "          -j INTEGER  : Number of threads to use, default 1" << std::endl;

        return -1;
    }

    // TODO:give better options
    // Lets handle all options
    int c;
    while((c = getopt(argc, argv, ":d:j:N:itoABCDE")) != -1){
        switch(c){
            case 'd':
                n_split = atoi(optarg);
                break;
            case 'j':
                n_threads = atoi(optarg);
            	break;
	        case 'N':
                Nstep = atoi(optarg);
                break;
            case 'i':
                InterlacedIter = true;
                break;
            case 't':
                TwoSideIter = true;
                break;
            case 'o':
                CosmicLaserIter = true;
                break;
            case 'A':
                DBoundary = true;
                break;
            case 'B':
                EBoundary = true;
                break;
            case 'C':
                CorrMapFlag = true;
                break;
            case 'W':
                WeightAverage = true;
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



//#ifdef _OPENMP
//    omp_set_num_threads(n_threads);
//#endif


    // Define input files for 2-side iteration
    std::vector<std::string> InputFiles1;
    std::vector<std::string> InputFiles2;

    // Define input file for interlaced iteration
    std::vector<std::string> InputFiles;

    unsigned int n_files = 0;

    for (int i = optind; i < argc; i++) {
        std::string filename(argv[i]);
        // check if file exists
        std::ifstream f(filename.c_str());
        if (!f.good()) {
            throw std::runtime_error(std::string("file does not exist: ") + filename);
        }

        if(TwoSideIter){
            TChain *tree = new TChain("lasers");
            tree->Add(filename.c_str());
            int side;
            tree->SetBranchAddress("side", &side);
            tree->Draw("side>>hside", "");
            TH1F *hside = (TH1F *) gDirectory->Get("hside");
            int LCS = hside->GetMean();
            delete tree;

            if (LCS == 1) {InputFiles1.push_back(filename); }
            else if(LCS==2){InputFiles2.push_back(filename); }
            else{ std::cerr << "The laser system is not labeled correctly." << std::endl; }
        }
        if(CosmicLaserIter){
            // The following is very dangerous! if there is "laser" and "cosmic" in the filename (path to the file)
            // The current version accept the file path contains "laser" but not "cosmic"
            std::string slaser ("aser");
            std::string scosmic ("osmic");
            if(filename.find(scosmic) != std::string::npos){
                std::cout<<"cosmic: "<<filename<<std::endl;
                InputFiles2.push_back(filename);
            }
            else if(filename.find(slaser) != std::string::npos) {
                std::cout<<"laser: "<<filename<<std::endl;
                InputFiles1.push_back(filename);
            }
            else{ std::cerr << "Laser or Cosmic? Check the file name." << std::endl; }
        }
        if(InterlacedIter){
            InputFiles.push_back(filename);
        }
    }

    if (TwoSideIter || CosmicLaserIter) {
        if(InputFiles1.empty() || InputFiles2.empty()){
            std::cerr << "Please provide the laser input data from 2 sides." << std::endl;
        }
    }
    else{
        if(InputFiles.empty()){
            std::cerr << "Please provide laser input data." << std::endl;
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

    std::stringstream ss_outfile;
    std::stringstream ss_Einfile;
    std::stringstream ss_Eoutfile;
    std::stringstream ss_D_outtxt;
    std::stringstream ss_E_outtxt;

    float float_max = std::numeric_limits<float>::max();
    ThreeVector<float> Unknown = {float_max, float_max, float_max};

    // Set the name for Dmap
    if (CorrMapFlag) {
        ss_outfile << "RecoCorr-N" << Nstep << "-S" << n_split << ".root";
        ss_D_outtxt << "RecoCorr-N" << Nstep << "-S" << n_split << ".txt";
    }
    if (!CorrMapFlag) {
        ss_outfile << "TrueDist-N" << Nstep << "-S" << n_split << ".root";
        ss_D_outtxt << "no_calib_usage_TrueDist-N" << Nstep << "-S" << n_split << ".txt";
    }

    // Name the input and output file name of E field calculation
    int NEinfile = 0;
    std::string Einfile;
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (".")) != NULL) {
        while ((ent = readdir (dir)) != NULL) {
            if(std::string(ent->d_name).compare(0,8,"RecoCorr")==0){
                NEinfile++;
                Einfile.assign(std::string(ent->d_name));
                ss_Einfile <<Einfile;
                ss_Eoutfile << "Emap-"<<Einfile.substr(9,Einfile.find_last_of(Einfile));
                ss_E_outtxt << "Emap"<< ".txt";
            }
        }
        closedir (dir);
    } else {
        std::cerr << "Local directory is not accessible." << std::endl;
    }

    if(NEinfile==1){
        std::ifstream ifile(ss_Einfile.str().c_str());
        if(ifile){
            std::cout<<"E field input file exists."<<std::endl;
        } else{
            std::cerr << "Please make sure there is one and only one 'RecoCorr*.root' file for E field calculation." << std::endl;
        }
    } else{
        std::cerr << "Please make sure there is one and only one 'RecoCorr*.root' file for E field calculation." << std::endl;
    }


    if (DoCorr) {
        std::vector<std::vector<ThreeVector<float>>> DisplMapsHolder;

        float float_max = std::numeric_limits<float>::max();
        ThreeVector<float> Empty = {float_max, float_max, float_max};

        // Read data and store it to a Laser object
        std::cout << "Reading data..." << std::endl;

        Laser TracksSample1;
        Laser TracksSample2;

        std::vector<Laser> LaserSets1;
        std::vector<Laser> LaserSets2;


        if(TwoSideIter || CosmicLaserIter){
            // Read the laser data from file to 'Laser'. The laser data from 2 sides are sperated
            TracksSample1 = ReadRecoTracks(InputFiles1);
            TracksSample2 = ReadRecoTracks(InputFiles2);
            // Split laser iteration samples into subsamples
            LaserSets1 = SplitTrackSet(TracksSample1, n_split);
            LaserSets2 = SplitTrackSet(TracksSample2, n_split);
        }

        else{
            // Read the laser data from file to 'Laser'.
            Laser FullTracks = ReadRecoTracks(InputFiles);
            // Seperate laser iteration sample
            TracksSample1 = InterlacedIterTrackSamples(FullTracks)[0];
            TracksSample2 = InterlacedIterTrackSamples(FullTracks)[1];
            // Split laser iteration samples into subsamples
            LaserSets1 = SplitTrackSet(TracksSample1, n_split);
            LaserSets2 = SplitTrackSet(TracksSample2, n_split);
        }

        int Mapsize = DetectorResolution[0]*DetectorResolution[1]*DetectorResolution[2];
//        std::pair<ThreeVector<float >, ThreeVector<float>>
//                PairIni = std::make_pair(ThreeVector<float>(0., 0., 0.),ThreeVector<float>(0., 0., 0.));
//
//        std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>> DisplacementMap(Mapsize, PairIni);
        std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>> DisplacementMap;
        DisplacementMap.reserve(Mapsize);

//        std::vector<ThreeVector<float>> DisplacementMap(DisplMapsHolder.front().size(),
//                                                        ThreeVector<float>(0., 0., 0.));
        if(WeightAverage){

            // Calculate track displacement
//            std::pair<Laser, Laser> LaserWithDisp = DispLaserIteration(Nstep, TracksSample1, TracksSample2, CorrMapFlag);
            std::pair<Laser, Laser> LaserWithDisp = DispLaserIteration(Nstep, LaserSets1[0], LaserSets2[0], CorrMapFlag);

            // TODO: Let's hope merge function is alright!
            Laser LaserCorrected = MergeLaser(LaserWithDisp.first, LaserWithDisp.second);

            std::cout << "Meshing for weighted mean in voxels" <<  std::endl;

            auto MeshforGrid = MeshVoxel(LaserCorrected.GetTrackSet(),Detector);

            std::cout << "Calculate weighted mean in voxels" <<  std::endl;

            DisplacementMap = AveragebyDistance(MeshforGrid, Detector);

        }

        if(!WeightAverage) {

            // Now we loop over each individual set and compute the displacement vectors.
            // TODO: This could be parallelized

//#pragma omp parallel for

            for (unsigned int set = 0; set < n_split; set++) {

                // The disadvantage is the LaserRecoOrigin will be discard after the calculation of this set
                Laser LaserRecoOrigin1 = LaserSets1[set];
                Laser LaserRecoOrigin2 = LaserSets2[set];

                std::cout << "Processing subset " << set << "/" << n_split << "... " << std::endl;

                // Calculate track displacement
                std::cout << " [" << set << "] Find track displacements... " << std::endl;

                std::pair<Laser, Laser> LaserWithDisp = DispLaserIteration(Nstep, LaserSets1[set], LaserSets2[set],
                                                                           CorrMapFlag);

                std::cout << "Time after N-step correction" << std::difftime(std::time(NULL), timer) << " s"
                          << std::endl;

//            // Merge 2 Laser samples with displacement vector for mesh and iteration
//            Laser LaserRecoOrigin = MergeLaser(LaserRecoOrigin1, LaserRecoOrigin2);
//            Laser LaserCorrected = MergeLaser(LaserWithDisp.first, LaserWithDisp.second);

                //Add anode information (no distortion) into Laser track sets
//            LaserRecoOrigin.AppendTrack(Anode(Detector));
//            LaserCorrected.AppendTrack(Anode(Detector));

                if (DBoundary) {
                    LaserRecoOrigin1.AppendTrack(Anode(Detector));
                    LaserRecoOrigin2.AppendTrack(Anode(Detector));
                    LaserWithDisp.first.AppendTrack(Anode(Detector));
                    LaserWithDisp.second.AppendTrack(Anode(Detector));
                }


                std::cout << " [" << set << "] Generate mesh..." << std::endl;

//            Delaunay MeshMap;
                Delaunay MeshMap1;
                Delaunay MeshMap2;

                std::cout << "Time after mesh " << std::difftime(std::time(NULL), timer) << " s" << std::endl;

                // The correction map is built on the mesh of reconstructed position which is the origin LaserSets
                if (CorrMapFlag) {
                    MeshMap1 = TrackMesher(LaserRecoOrigin1.GetTrackSet());
                    MeshMap2 = TrackMesher(LaserRecoOrigin2.GetTrackSet());
                }

                    // The distortion map is built on the mesh of true position which is moved LaserSets
                else {
                    MeshMap1 = TrackMesher(LaserWithDisp.first.GetTrackSet());
                    MeshMap2 = TrackMesher(LaserWithDisp.second.GetTrackSet());
                }

                // Interpolate Displacement Map (regularly spaced grid)
                std::cout << "Start interpolation..." << std::endl;
                // LaserSets are now sitting on the true position, LaserRecoOrigin are sitting on the reco position

                // The correction map is based on reco space coord
                if (CorrMapFlag) {
                    DisplMapsHolder.push_back(
                            InterpolateMap(LaserWithDisp.first.GetTrackSet(), LaserRecoOrigin1.GetTrackSet(), MeshMap1,
                                           Detector));
                    DisplMapsHolder.push_back(
                            InterpolateMap(LaserWithDisp.second.GetTrackSet(), LaserRecoOrigin2.GetTrackSet(), MeshMap2,
                                           Detector));
                }

                    // The distortion map is based on true space coord
                else {
                    DisplMapsHolder.push_back(
                            InterpolateMap(LaserWithDisp.first.GetTrackSet(), LaserWithDisp.first.GetTrackSet(),
                                           MeshMap1, Detector));
                    DisplMapsHolder.push_back(
                            InterpolateMap(LaserWithDisp.second.GetTrackSet(), LaserWithDisp.second.GetTrackSet(),
                                           MeshMap2, Detector));
                }
            }

            // Now we go on to create an unified displacement map
//            std::vector<ThreeVector<float>> DisplacementMap(DisplMapsHolder.front().size(),
//                                                            ThreeVector<float>(0., 0., 0.));
            std::vector<float> Nvalid(DisplMapsHolder.front().size(), 0.);

            for (auto &SubMap: DisplMapsHolder) {
                for (unsigned int idx = 0; idx < DisplacementMap.size(); idx++) {
                    if (SubMap[idx] != Empty) {
                        DisplacementMap[idx].first = DisplacementMap[idx].first + SubMap[idx];
                        Nvalid[idx]++;
                    }
                }
            }

            for (unsigned int idx = 0; idx < DisplacementMap.size(); idx++) {
                if (Nvalid[idx] == 0) {
                    // Set those bin with non valid number into float max again
                    DisplacementMap[idx].first = {float_max, float_max, float_max};
                } else {
                    DisplacementMap[idx].first = DisplacementMap[idx].first / Nvalid[idx];
                }
            }


//            // Loop the displacement map in the form of vector
//            for (int idx = 0; idx < DisplacementMap.size(); idx++){
//
//                std::vector<ThreeVector<float>> BinStatistics;
//                ThreeVector<float> BinAverage;
//                ThreeVector<float> BinStd;
//                int Nvalid = 0;
//
//                // Loop the bins in submaps to calculate the bin averaging displacement
//                for (int NsubMap = 0; NsubMap < DisplMapsHolder.size(); NsubMap++){
//                    if(DisplMapsHolder[NsubMap][idx] != Empty){
//                        // If a submap bin is not empty, push back to a profile vector
//                        BinStatistics.push_back(DisplMapsHolder[NsubMap][idx]);
//                        BinAverage += BinStatistics.back();
//                        Nvalid++;
//                    }
//                }
//                // Averaging displacement in a bin
//                BinAverage = BinAverage / (float) Nvalid;
//
//                // Loop the profile vector to calculate the standard deviation
//                for (int binValid = 0; binValid < Nvalid; binValid++){
//                    for(int k = 0; k < 3; k++){
//                        BinStd[k] += (BinStatistics[binValid][k]-BinAverage[k])* (BinStatistics[binValid][k]-BinAverage[k]);
//                    }
//                }
//                // Standard deviation of the displacement in a bin
//                for(int k = 0; k < 3; k++){
//                    BinStd[k] = sqrtf(BinStd[k] / (float) Nvalid);
//                }
//
//                // Construct 1d displacement map
//                std::pair<ThreeVector<float >, ThreeVector<float>> BinInfo = std::make_pair(BinAverage,BinStd);
//                DisplacementMap.push_back(BinInfo);
//
//            }


//            // Fill displacement map into TH3 histograms and write them to file
//            std::cout << "Write to File ..." << std::endl;
//            WriteRootFile(DisplacementMap, Detector, ss_outfile.str());
        }

        // Fill displacement map into TH3 histograms and write them to root and txt file
        std::cout << "Write to File ..." << std::endl;
        WriteRootFileDeviation(DisplacementMap, Detector, ss_outfile.str());
        WriteTextFile(DisplacementMap,ss_D_outtxt.str());
    }

    // The Emap calculation works when the input is correction map
    if (DoEmap) {

        // The vector of Position and En must have the exactly the same index to make the interpolation (EInterpolateMap()) work
        if(EBoundary){

            auto EfieldXYZ = EfieldXYZwithBoundary(Detector, cryoTemp, E0, v0, ss_Einfile.str().c_str());
            std::vector<float> Ex = std::get<0>(EfieldXYZ);
            std::vector<float> Ey = std::get<1>(EfieldXYZ);
            std::vector<float> Ez = std::get<2>(EfieldXYZ);
            std::vector<ThreeVector<float>> PositionX = std::get<3>(EfieldXYZ);
            std::vector<ThreeVector<float>> PositionY = std::get<3>(EfieldXYZ);
            std::vector<ThreeVector<float>> PositionZ = std::get<3>(EfieldXYZ);

            // Create mesh for Emap
            std::cout << "Generate mesh for E field..." << std::endl;
            xDelaunay EMeshX = Mesher(PositionX, Detector);
            xDelaunay EMeshY = Mesher(PositionY, Detector);
            xDelaunay EMeshZ = Mesher(PositionZ, Detector);

            // Interpolate E Map (regularly spaced grid)
            std::cout << "Start interpolation the E field..." << std::endl;
            std::vector<ThreeVector<float>> EMapXYZ = EcompInterpolateMap(Ex, PositionX, EMeshX,
                                                                          Ey, PositionY, EMeshY,
                                                                          Ez, PositionZ, EMeshZ,
                                                                          Detector, EMapResolution);

            // Fill displacement map into TH3 histograms and write them to file
            std::cout << "Write Emap to File ..." << std::endl;
            WriteEmapRoot(EMapXYZ, Detector, EMapResolution, E0, ss_Eoutfile.str());

        }
        else{
            // The vector of Position and En must have the exactly the same index to make the interpolation (EInterpolateMap()) work
            auto E_field = Efield(Detector, cryoTemp, E0, v0, ss_Einfile.str().c_str());
            std::vector<ThreeVector<float>> En = E_field.first;
            std::vector<ThreeVector<float>> Position = E_field.second;

            // Create mesh for Emap
            std::cout << "Generate mesh for E field..." << std::endl;
            xDelaunay EMesh = Mesher(Position, Detector);

            // Interpolate E Map (regularly spaced grid)
            std::cout << "Start interpolation the E field..." << std::endl;
            std::vector<ThreeVector<float>> EMap = EInterpolateMap(En, Position, EMesh, Detector, EMapResolution);

            // Fill displacement map into TH3 histograms and write them to file
            std::cout << "Write Emap to File ..." << std::endl;
            WriteEmapRoot(EMap, Detector, EMapResolution, E0, ss_Eoutfile.str());
        }



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

            /*
            // This here sorts the tracks by their distance to the EntryPoint. The algorithm uses a lambda
            // It will compare the distance to the EntryPoint of two vector entries A & B
            std::sort(TrackSamples.begin(), TrackSamples.end(), [&EntryPoint](TVector3 A, TVector3 B) {
                          A -= EntryPoint;
                          B -= EntryPoint;
                          // Here only the squared distance was used to avoid costly sqrt operations
                          return A.Mag2() > B.Mag2();
                      }
            );
            */

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
                   std::string OutputFilename) {
    // Store TPC properties which are important for the TH3 generation

    ThreeVector<unsigned long> Resolution = TPCVolume.GetDetectorResolution();
    ThreeVector<float> MinimumCoord = TPCVolume.GetMapMinimum();
    ThreeVector<float> MaximumCoord = TPCVolume.GetMapMaximum();
    ThreeVector<float> Unit = {TPCVolume.GetDetectorSize()[0] / (Resolution[0] - 1),
                               TPCVolume.GetDetectorSize()[1] / (Resolution[1] - 1),
                               TPCVolume.GetDetectorSize()[2] / (Resolution[2] - 1)};

    // Initialize all TH3F
    std::vector<TH3F> RecoDisplacement;
    RecoDisplacement.push_back(TH3F("Reco_Displacement_X", "Reco Displacement X",
                                    Resolution[0], MinimumCoord[0] - Unit[0] * 0.5, MaximumCoord[0] + Unit[0] * 0.5,
                                    Resolution[1], MinimumCoord[1] - Unit[1] * 0.5, MaximumCoord[1] + Unit[1] * 0.5,
                                    Resolution[2], MinimumCoord[2] - Unit[2] * 0.5, MaximumCoord[2] + Unit[2] * 0.5));
    RecoDisplacement.push_back(TH3F("Reco_Displacement_Y", "Reco Displacement Y",
                                    Resolution[0], MinimumCoord[0] - Unit[0] * 0.5, MaximumCoord[0] + Unit[0] * 0.5,
                                    Resolution[1], MinimumCoord[1] - Unit[1] * 0.5, MaximumCoord[1] + Unit[1] * 0.5,
                                    Resolution[2], MinimumCoord[2] - Unit[2] * 0.5, MaximumCoord[2] + Unit[2] * 0.5));
    RecoDisplacement.push_back(TH3F("Reco_Displacement_Z", "Reco Displacement Z",
                                    Resolution[0], MinimumCoord[0] - Unit[0] * 0.5, MaximumCoord[0] + Unit[0] * 0.5,
                                    Resolution[1], MinimumCoord[1] - Unit[1] * 0.5, MaximumCoord[1] + Unit[1] * 0.5,
                                    Resolution[2], MinimumCoord[2] - Unit[2] * 0.5, MaximumCoord[2] + Unit[2] * 0.5));

    // Loop over all xbins
    for (unsigned xbin = 0; xbin < Resolution[0] ; xbin++) {
        // Loop over all ybins
        for (unsigned ybin = 0; ybin < Resolution[1] ; ybin++) {
            // Loop over all zbins
            for (unsigned zbin = 0; zbin < Resolution[2] ; zbin++) {
                // Loop over all coordinates

                for (unsigned coord = 0; coord < 3; coord++) {
                    // Fill interpolated grid points into histograms
                    RecoDisplacement[coord].SetBinContent(xbin + 1, ybin + 1, zbin + 1,
                                                          InterpolationData[zbin + (Resolution[2] * (ybin + Resolution[1] * xbin))][coord]);
                    // Remember, the range of the hist bin is (1, nbins), while when we fill the vector, it starts from 0. (0,nbins-1)
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

void WriteRootFileDeviation(std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>> &InterpolationData,
                            TPCVolumeHandler &TPCVolume, std::string OutputFilename) {
    // Store TPC properties which are important for the TH3 generation

    ThreeVector<unsigned long> Resolution = TPCVolume.GetDetectorResolution();
    ThreeVector<float> MinimumCoord = TPCVolume.GetMapMinimum();
    ThreeVector<float> MaximumCoord = TPCVolume.GetMapMaximum();
    ThreeVector<float> Unit = {TPCVolume.GetDetectorSize()[0] / (Resolution[0] - 1),
                               TPCVolume.GetDetectorSize()[1] / (Resolution[1] - 1),
                               TPCVolume.GetDetectorSize()[2] / (Resolution[2] - 1)};

    // Initialize all TH3F
    std::vector<TH3F> RecoDisplacement(6, TH3F("Reco_Displacement", "Reco Displacement",
                                    Resolution[0], MinimumCoord[0] - Unit[0] * 0.5, MaximumCoord[0] + Unit[0] * 0.5,
                                    Resolution[1], MinimumCoord[1] - Unit[1] * 0.5, MaximumCoord[1] + Unit[1] * 0.5,
                                    Resolution[2], MinimumCoord[2] - Unit[2] * 0.5, MaximumCoord[2] + Unit[2] * 0.5));

    RecoDisplacement[0].SetNameTitle("Reco_Displacement_X", "Reco Displacement X");
    RecoDisplacement[1].SetNameTitle("Reco_Displacement_Y", "Reco Displacement Y");
    RecoDisplacement[2].SetNameTitle("Reco_Displacement_Z", "Reco Displacement Z");
    RecoDisplacement[3].SetNameTitle("Reco_Displacement_X_Error", "Reco Deviation of Displacement X");
    RecoDisplacement[4].SetNameTitle("Reco_Displacement_Y_Error", "Reco Deviation of Displacement X");
    RecoDisplacement[5].SetNameTitle("Reco_Displacement_Z_Error", "Reco Deviation of Displacement X");

    // Loop over all xbins
    for (unsigned xbin = 0; xbin < Resolution[0] ; xbin++) {
        // Loop over all ybins
        for (unsigned ybin = 0; ybin < Resolution[1] ; ybin++) {
            // Loop over all zbins
            for (unsigned zbin = 0; zbin < Resolution[2] ; zbin++) {
                // Loop over all coordinates

                for (unsigned coord = 0; coord < 3; coord++) {
                    // Fill interpolated grid points into histograms
                    // the range of the hist bin is (1, nbins), while when we fill the vector, it starts from 0. (0,nbins-1)
                    RecoDisplacement[coord].SetBinContent(xbin + 1, ybin + 1, zbin + 1,
                                                          InterpolationData[zbin + (Resolution[2] * (ybin + Resolution[1] * xbin))].first[coord]);
                    RecoDisplacement[coord+3].SetBinContent(xbin + 1, ybin + 1, zbin + 1,
                                                          InterpolationData[zbin + (Resolution[2] * (ybin + Resolution[1] * xbin))].second[coord]);
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


void WriteTextFile(std::vector<std::pair<ThreeVector<float >, ThreeVector<float>>> &InterpolationData,
                   std::string OutputFilename) {
    // Initialize stream to file
    std::ofstream OutputFile;

    // Open output file
    OutputFile.open(OutputFilename.c_str(), std::ios::out);

    // Loop over all interpolated data points
    for (unsigned entry = 0; entry < InterpolationData.size(); entry++) {
        // Write every point into a seperate line
        OutputFile << entry <<"\t"
                   << InterpolationData[entry].first[0] <<"\t" << InterpolationData[entry].first[1] <<"\t"
                   << InterpolationData[entry].first[2] <<"\t"
                   << InterpolationData[entry].second[0] <<"\t" << InterpolationData[entry].second[1] <<"\t"
                   << InterpolationData[entry].second[2];
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

// Write Emap into TH3 and store in root file
void WriteEmapRoot(std::vector<ThreeVector<float>> &Efield, TPCVolumeHandler &TPCVolume,
                   ThreeVector<unsigned long> Resolution, float E0, std::string OutputFilename) {
    // Store TPC properties which are important for the TH3 generation
    float float_max = std::numeric_limits<float>::max();
    ThreeVector<float> MinimumCoord = TPCVolume.GetMapMinimum();
    ThreeVector<float> MaximumCoord = TPCVolume.GetMapMaximum();
    ThreeVector<float> Unit = {TPCVolume.GetDetectorSize()[0] / (Resolution[0] - 1),
                               TPCVolume.GetDetectorSize()[1] / (Resolution[1] - 1),
                               TPCVolume.GetDetectorSize()[2] / (Resolution[2] - 1)};

    // Initialize all TH3F
    std::vector<TH3F> Emap;
    Emap.push_back(TH3F("Emap_X", "E field map X",
                        Resolution[0], MinimumCoord[0] - Unit[0] * 0.5, MaximumCoord[0] + Unit[0] * 0.5,
                        Resolution[1], MinimumCoord[1] - Unit[1] * 0.5, MaximumCoord[1] + Unit[1] * 0.5,
                        Resolution[2], MinimumCoord[2] - Unit[2] * 0.5, MaximumCoord[2] + Unit[2] * 0.5));
    Emap.push_back(TH3F("Emap_Y", "E field map Y",
                        Resolution[0], MinimumCoord[0] - Unit[0] * 0.5, MaximumCoord[0] + Unit[0] * 0.5,
                        Resolution[1], MinimumCoord[1] - Unit[1] * 0.5, MaximumCoord[1] + Unit[1] * 0.5,
                        Resolution[2], MinimumCoord[2] - Unit[2] * 0.5, MaximumCoord[2] + Unit[2] * 0.5));
    Emap.push_back(TH3F("Emap_Z", "E field map Z",
                        Resolution[0], MinimumCoord[0] - Unit[0] * 0.5, MaximumCoord[0] + Unit[0] * 0.5,
                        Resolution[1], MinimumCoord[1] - Unit[1] * 0.5, MaximumCoord[1] + Unit[1] * 0.5,
                        Resolution[2], MinimumCoord[2] - Unit[2] * 0.5, MaximumCoord[2] + Unit[2] * 0.5));


    // the loop should be consistent to the one in the EInterpolateMap()
    for (unsigned xbin = 0; xbin < Resolution[0]; xbin++) {
        for (unsigned ybin = 0; ybin < Resolution[1]; ybin++) {
            for (unsigned zbin = 0; zbin < Resolution[2]; zbin++) {
//                if((Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][0] > 0.5*float_max ||
//                   Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][1] > 0.5*float_max ||
//                   Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][2] > 0.5*float_max) &&
//                   xbin == 0)
//                {
////                    std::cout<<"x: "<<xbin<<"; y: "<<ybin<<"; z: "<<zbin<<" !Flag!"<<std::endl;
//                    ThreeVector<unsigned long> Coord = {xbin, ybin, zbin};
//                    Emap[0].SetBinContent(xbin + 1, ybin + 1, zbin + 1, EdgeEx(Efield, Resolution, Unit, E0, Coord));
//                    Emap[1].SetBinContent(xbin + 1, ybin + 1, zbin + 1, 0);
//                    Emap[2].SetBinContent(xbin + 1, ybin + 1, zbin + 1, 0);
//                    std::cout<<"Maxwell...xbin: "<<xbin<<", ybin: "<<ybin<<", zbin: "<<zbin<<", Ex: "<<EdgeEx(Efield, Resolution, Unit, E0, Coord)<<std::endl;
//                }
//                else if(xbin ==0){
//                    Emap[0].SetBinContent(xbin + 1, ybin + 1, zbin + 1, -999);
//                    Emap[1].SetBinContent(xbin + 1, ybin + 1, zbin + 1, -999);
//                    Emap[2].SetBinContent(xbin + 1, ybin + 1, zbin + 1, -999);
//                }
//                else{
//                    if(xbin == 0){
//                        std::cout<<"Mesh-Interpolation...xbin: "<<xbin<<", ybin: "<<ybin<<", zbin: "<<zbin
//                                 <<", Ex: "<<Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][0]
//                                 <<", Ey: "<<Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][1]
//                                 <<", Ez: "<<Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][2]
//                                 <<std::endl;
//                    }
                    // Loop over all coordinates dx,dy,dz
                    for (unsigned coord = 0; coord < 3; coord++) {
                        // Fill interpolated grid points into histograms. bin=0 is underflow, bin = nbin+1 is overflow
                        Emap[coord].SetBinContent(xbin + 1, ybin + 1, zbin + 1,
                                                  Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][coord]);
//                        if(Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][coord]>20){
//                            std::cout<<"x: "<<xbin +1 <<"; y: "<<ybin+1<<"; z: "<<zbin+1<<"Efield[ "<<coord<<"]: "
//                                     <<Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][coord]
//                                     <<std::endl;
//                        }
                    } // end coordinate loop
//                }
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
