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

// Initialize functions defined below
//bool Twolasersys(int argc, char** argv);
Laser ReadRecoTracks(int argc, char** argv);
void WriteRootFile(std::vector<ThreeVector<float>>& , TPCVolumeHandler&);
void WriteTextFile(std::vector<ThreeVector<float>>&);
void LaserInterpThread(Laser&, const Laser&, const Delaunay&);
std::vector<Laser> ReachedExitPoint(const Laser&, float);
std::vector<ThreeVector<float>> Elocal(TPCVolumeHandler& );
std::vector<ThreeVector<float>> Eposition(TPCVolumeHandler& );
void WriteEmapRoot(std::vector<ThreeVector<float>>& Efield, TPCVolumeHandler& TPCVolume);

// Set if the output displacement map is correction map (on reconstructed coordinate) or distortion map (on true coordinate)
// By default set it as correction map so we could continue calculate the E field map
bool CorrMapFlag = true;
bool DoCorr = true;
bool DoEmap = true;

// Main function
int main(int argc, char** argv) {
    // Start timer, just because it's nice to know how long shit takes
    time_t timer;
    std::time(&timer);

    // If there are to few input arguments abord
    if (argc < 2) {
        std::cerr << "ERROR: Too few arguments, use ./LaserCal <input file names>" << std::endl;
        return -1;
    }

    // Choose detector dimensions, coordinate system offset and resolutions
    ThreeVector<float> DetectorSize = {256.04, 232.5, 1036.8};
    ThreeVector<float> DetectorOffset = {0.0, -DetectorSize[1] / static_cast<float>(2.0), 0.0};
    ThreeVector<unsigned long> DetectorResolution = {26, 26, 101};
    // Create the detector volume
    TPCVolumeHandler Detector(DetectorSize, DetectorOffset, DetectorResolution);

    if(DoCorr){

        // Read data and store it to a Laser object
        std::cout << "Reading data..." << std::endl;
        Laser LaserTrackSet = ReadRecoTracks(argc, argv);

        // Calculate track displacement
        std::cout << "Find track displacements... " << std::endl;

        if (CorrMapFlag) {
            // Choose displacement algorithm (available so far: TrackDerivative, ClosestPoint, or LinearStretch)
            // Suggestion: stay with ClosestPoint Algorithm
            LaserTrackSet.CalcDisplacement(LaserTrack::ClosestPointCorr);

            // Now the laser data are based on the reconstructed coordinate. For CORRECTION MAP, they are good enough.
            // No need to prepare to set true coordinate
        }

        // Caculating now the displacement map of distortion based on the true coordinates
        // Remember to turn the "CorrMapFlag" off
        if (!CorrMapFlag) {
            // Choose displacement algorithm (available so far: TrackDerivative, ClosestPoint, or LinearStretch)
            // Suggestion: stay with ClosestPoint Algorithm
            LaserTrackSet.CalcDisplacement(LaserTrack::ClosestPointDist);

            // Now the laser tracks are based on the reconstructed coordinate. If require DISTORTION MAP as output, set the base on the true coordinate
            // Add displacement to reconstructed track to change to detector coordinates (only for map generation)
            LaserTrackSet.AddCorrectionToReco();
        }


        // Create delaunay mesh
        std::cout << "Generate mesh..." << std::endl;
        Delaunay Mesh = TrackMesher(LaserTrackSet.GetTrackSet());

        // Interpolate Displacement Map (regularly spaced grid)
        std::cout << "Start interpolation..." << std::endl;
        std::vector<ThreeVector<float>> DisplacementMap = InterpolateMap(LaserTrackSet.GetTrackSet(), Mesh, Detector);

        // Fill displacement map into TH3 histograms and write them to file
        std::cout << "Write to File ..." << std::endl;
        WriteRootFile(DisplacementMap, Detector);
    }

    // Start the session to calculate E map
    // The Emap calculation works when the input is correction map
    if(CorrMapFlag && DoEmap){
        // The vector of Position and En must have the exactly the same index to make the interpolation (EInterpolateMap()) work
        std::vector<ThreeVector<float>> Position = Eposition(Detector);
        std::vector<ThreeVector<float>> En = Elocal(Detector);

        // Create mesh for Emap
        std::cout << "Generate mesh for E field..." << std::endl;
        xDelaunay EMesh = Mesher(Position, Detector);

        // Interpolate E Map (regularly spaced grid)
        std::cout << "Start interpolation the E field..." << std::endl;
        std::vector<ThreeVector<float>> EMap = EInterpolateMap(En, Position, EMesh, Detector);

        // Fill displacement map into TH3 histograms and write them to file
        std::cout << "Write Emap to File ..." << std::endl;
        WriteEmapRoot(EMap,Detector);
    }

    std::cout << "End of program after "<< std::difftime(std::time(NULL),timer) << " s" << std::endl;


} // end main

// To check if the input root files contain the laser data from two laser system (two side)
//bool Twolasersys(int argc, char** argv)
//{
//    if(argc != 3){
//        std::cerr << "ERROR: Not right arguments. Please use ./FieldCal <name_LCS1.root> <name_LCS2.root>" << std::endl;
//        return false; //0
//    }
//    if(argc == 3){
//        int LCS1,LCS2;
//        TChain* tree1 = new TChain("laser1");
//        tree1->Add(argv[1]);
//        tree1->SetBranchAddress("LCS1",&LCS1);
//        TH1F* h1;
//        tree1->Draw("LCS1>>h1","");
//        LCS1 = h1->GetMean();
//
//        TChain* tree2 = new TChain("laser2");
//        tree2->Add(argv[2]);
//        tree2->SetBranchAddress("LCS2",&LCS2);
//        TH1F* h2;
//        tree2->Draw("LCS1>>h2","");
//        LCS2 = h2->GetMean();
//
//        // A rough check here. Risk that one single input root file has mixed data from two laser system
//        if((LCS1-1.5)*(LCS2-1.5)<0){ return true; }
//        else{return false;}
//    }
//}

Laser ReadRecoTracks(int argc, char** argv)
{
    // Create Laser (collection of laser tracks) this will be the returned object
    Laser TrackSelection;
    
    // Initialize read variables, the pointers for more complex data structures 
    // are very important for Root. Rene Brun in hell (do you see what I did there?)
    int EventNumber;
    
    std::vector<TVector3> TrackSamples;
    std::vector<TVector3>* pTrackSamples = &TrackSamples;
    
    TVector3 EntryPoint;
    TVector3* pEntryPoint = &EntryPoint;
    TVector3 ExitPoint;
    TVector3* pExitPoint =&ExitPoint;
    
    // Open TChains to store all trees
    TChain* LaserInfoTree = new TChain("lasers");
    TChain* RecoTrackTree = new TChain("tracks");
    
    // Loop through all input files and add them to the TChain
    for(int arg = 1; arg < argc; arg++)
    {
        // Open input file and add to TChains
        LaserInfoTree->Add(argv[arg]);
        RecoTrackTree->Add(argv[arg]);
    }
    
    // Assigne branch addresses
    LaserInfoTree->SetBranchAddress("entry",&pEntryPoint);
    LaserInfoTree->SetBranchAddress("exit", &pExitPoint);
    RecoTrackTree->SetBranchAddress("track",&pTrackSamples);
    RecoTrackTree->SetBranchAddress("event", &EventNumber);
    
    // Only start read out when both trees have the same amount of entries 
    if(LaserInfoTree->GetEntries() == RecoTrackTree->GetEntries())
    {
        // Loop over all tree entries
        for(Size_t tree_index = 0; tree_index < RecoTrackTree->GetEntries(); tree_index++)
        {
            // Get tree entries of both trees
            LaserInfoTree->GetEntry(tree_index);
            RecoTrackTree->GetEntry(tree_index);
            

            // This here sorts the tracks by their distance to the EntryPoint. The algorithm uses a lambda
            // It will compare the distance to the EntryPoint of two vector entries A & B
            std::sort(TrackSamples.begin(),TrackSamples.end(), [&EntryPoint](TVector3 A, TVector3 B)
                {
                    A -= EntryPoint;
                    B -= EntryPoint;
                    // Her only the squared distance was used to avoid costly sqrt operations
                    return A.Mag2() > B.Mag2();
                }    
            );
            
            // This step will erase all double entries. First std::unique shifts every double to the end
            // of the vector and gives back the new end point of the data set. After that we erase the HistRange
            // between this new end and the real end of the vector
            TrackSamples.erase( std::unique(TrackSamples.begin(),TrackSamples.end()), TrackSamples.end());
            
            // Add new track to Laser TrackSelection
            TrackSelection.AppendTrack(LaserTrack(EntryPoint,ExitPoint,TrackSamples));
        }
    }
    else // If the trees don't have the same amount of entries, through error (I know not propper error handling)
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


void WriteRootFile(std::vector<ThreeVector<float>>& InterpolationData, TPCVolumeHandler& TPCVolume)
{ 
    // Store TPC properties which are important for the TH3 generation
    ThreeVector<unsigned long> Resolution = TPCVolume.GetDetectorResolution();
    ThreeVector<float> MinimumCoord = TPCVolume.GetMapMinimum();
    ThreeVector<float> MaximumCoord = TPCVolume.GetMapMaximum();
  
    // Initialize all TH3F
    std::vector<TH3F> RecoDisplacement;
    RecoDisplacement.push_back(TH3F("Reco_Displacement_X","Reco Displacement X",Resolution[0],MinimumCoord[0],MaximumCoord[0],Resolution[1],MinimumCoord[1],MaximumCoord[1],Resolution[2],MinimumCoord[2],MaximumCoord[2]));
    RecoDisplacement.push_back(TH3F("Reco_Displacement_Y","Reco Displacement Y",Resolution[0],MinimumCoord[0],MaximumCoord[0],Resolution[1],MinimumCoord[1],MaximumCoord[1],Resolution[2],MinimumCoord[2],MaximumCoord[2]));
    RecoDisplacement.push_back(TH3F("Reco_Displacement_Z","Reco Displacement Z",Resolution[0],MinimumCoord[0],MaximumCoord[0],Resolution[1],MinimumCoord[1],MaximumCoord[1],Resolution[2],MinimumCoord[2],MaximumCoord[2]));
  

    // Loop over all xbins
    for(unsigned xbin = 0; xbin < TPCVolume.GetDetectorResolution()[0]; xbin++)
    {
        // Loop over all ybins
        for(unsigned ybin = 0; ybin < TPCVolume.GetDetectorResolution()[1]; ybin++) 
        {
            // Loop over all zbins
            for(unsigned zbin = 0; zbin < TPCVolume.GetDetectorResolution()[2]; zbin++)
            {
                // Loop over all coordinates
                for(unsigned coord = 0; coord < 3; coord++)
                {
                    // Fill interpolated grid points into histograms
                    // bin=0 is underflow, bin = nbin+1 is overflow
                    RecoDisplacement[coord].SetBinContent(xbin+1,ybin+1,zbin+1, InterpolationData[zbin+( TPCVolume.GetDetectorResolution()[2]*(ybin+TPCVolume.GetDetectorResolution()[1]*xbin) )][coord]);
                    // It's equivalent to the following expression
                    // Remember, the range of the hist bin is (1, nbins), while when we fill the vector, it starts from 0. (0,nbins-1)
                    // RecoDisplacement[coord].SetBinContent(xbin+1,ybin+1,zbin+1, InterpolationData[zbin+ybin*TPCVolume.GetDetectorResolution()[2]+xbin*TPCVolume.GetDetectorResolution()[1]*TPCVolume.GetDetectorResolution()[2]][coord]);
                } // end coordinate loop
            } // end zbin loop
        } // end ybin loop
    } // end zbin loop
    
    // Open and recreate output file
    TFile OutputFile("RecoCorr.root", "recreate");
    
    // Loop over space coordinates
    for(unsigned coord = 0; coord < RecoDisplacement.size(); coord++)
    {
        // Write every TH3 map into file
        RecoDisplacement[coord].Write();
    }
    
    // Close output file and clean up
    OutputFile.Close();
    gDirectory->GetList()->Delete();
}


void WriteTextFile(std::vector<ThreeVector<float>>& InterpolationData)
{
    // Initialize stream to file
    std::ofstream OutputFile;
  
    // Open output file
    OutputFile.open("Reco.txt", std::ios::out);
  
    // Loop over all interpolated data points
    for(unsigned entry = 0; entry < InterpolationData.size(); entry++)
    {
        // Write every point into a seperate line
        OutputFile << InterpolationData[entry][0] << InterpolationData[entry][1] << InterpolationData[entry][2];
    }

    // Close file
    OutputFile.close();
} // WriteRootFile


// This is the multi-threading interpolation function. Just hangs out here, for legacy purposes 
void LaserInterpThread(Laser& LaserTrackSet, const Laser& InterpolationLaser, const Delaunay& InterpolationMesh)
{
  LaserTrackSet.InterpolateTrackSet(InterpolationLaser, InterpolationMesh);
} // LaserInterpThread

// Split the laser track set into tracks that reached the expected exit point (within a configurable region) and others.
// First entry of the return vector is tracks that reach the exit point, second is the ones that do not reach it.


// The root file does not have to be the argument
std::vector<ThreeVector<float>> Elocal(TPCVolumeHandler& TPCVolume)
{
    TFile *InFile = new TFile("RecoCorr.root","READ");

    TH3F *Dx = (TH3F*) InFile->Get("Reco_Displacement_X");
    TH3F *Dy = (TH3F*) InFile->Get("Reco_Displacement_Y");
    TH3F *Dz = (TH3F*) InFile->Get("Reco_Displacement_Z");

    float DetectorReso[3]={TPCVolume.GetDetectorSize()[0] /TPCVolume.GetDetectorResolution()[0],TPCVolume.GetDetectorSize()[1]/TPCVolume.GetDetectorResolution()[1],TPCVolume.GetDetectorSize()[2]/TPCVolume.GetDetectorResolution()[2]};
    float Delta_x = DetectorReso[0]; //cm
    float Ex = 273; // kV/cm

    std::vector<ThreeVector< float>> En(TPCVolume.GetDetectorResolution()[2] * TPCVolume.GetDetectorResolution()[1] * (TPCVolume.GetDetectorResolution()[0]-1));


    for(unsigned zbin = 0; zbin < TPCVolume.GetDetectorResolution()[2]; zbin++)
    {
        for(unsigned ybin = 0; ybin < TPCVolume.GetDetectorResolution()[1]; ybin++)
        {
            // the number of x bin in Emap is one less than in the displacement map
            // because we only consider the gap here
            for(unsigned xbin = 0; xbin < (TPCVolume.GetDetectorResolution()[0]-1); xbin++)
            {
                //xbin =1, x=5 close to the anode; xbin = Nx, x=250 close to the cathode
                // the corner has the weight of the bin. Do not move the weight to the geometry center!
                ThreeVector<float> RecoGrid(xbin * DetectorReso[0],ybin*DetectorReso[1]+TPCVolume.GetDetectorOffset()[1],zbin * DetectorReso[2]);
                ThreeVector<float> Dxyz = {(float) Dx->GetBinContent(xbin+1,ybin+1,zbin+1), (float) Dy->GetBinContent(xbin+1,ybin+1,zbin+1), (float) Dz->GetBinContent(xbin+1,ybin+1,zbin+1)};
                ThreeVector<float> True = RecoGrid + Dxyz;

                ThreeVector<float> RecoGrid_next((xbin+1)*DetectorReso[0],ybin*DetectorReso[1]+TPCVolume.GetDetectorOffset()[1],zbin*DetectorReso[2]);
                ThreeVector<float> Dxyz_next = {(float) Dx->GetBinContent(xbin+2,ybin+1,zbin+1),(float) Dy->GetBinContent(xbin+2,ybin+1,zbin+1),(float) Dz->GetBinContent(xbin+2,ybin+1,zbin+1)};
                ThreeVector<float> True_next = RecoGrid_next + Dxyz_next;

                ThreeVector<float> Rn = True_next - True;
                En[xbin+ybin*(TPCVolume.GetDetectorResolution()[0]-1)+zbin*(TPCVolume.GetDetectorResolution()[0]-1)*TPCVolume.GetDetectorResolution()[1]] = Ex / Delta_x * Rn;
            }
        }
    }
    return En;
}

std::vector<ThreeVector<float>> Eposition(TPCVolumeHandler& TPCVolume)
{
    ThreeVector<unsigned long> Resolution = TPCVolume.GetDetectorResolution();

    TFile *InFile = new TFile("RecoCorr.root","READ");

    TH3F *Dx = (TH3F*) InFile->Get("Reco_Displacement_X");
    TH3F *Dy = (TH3F*) InFile->Get("Reco_Displacement_Y");
    TH3F *Dz = (TH3F*) InFile->Get("Reco_Displacement_Z");

    float DetectorReso[3]={TPCVolume.GetDetectorSize()[0] /Resolution[0],TPCVolume.GetDetectorSize()[1]/Resolution[1],TPCVolume.GetDetectorSize()[2]/Resolution[2]};

    std::vector<ThreeVector<float>> Position(Resolution[2]*Resolution[1]*(Resolution[0]-1));

    // the position should be consistent to the one in the EInterpolateMap()
    for(unsigned zbin = 0; zbin < Resolution[2]; zbin++)
    {
        for(unsigned ybin = 0; ybin < Resolution[1]; ybin++)
        {
            // since we calculate Elocal by the gap, the number of x bin in Emap is one less than in the displacement map
            for(unsigned xbin = 0; xbin < (Resolution[0]-1); xbin++)
            {
                //xbin =1, x=5 close to the anode; xbin = Nx, x=250 close to the cathode
                ThreeVector<float> RecoGrid(xbin * DetectorReso[0],ybin*DetectorReso[1]+TPCVolume.GetDetectorOffset()[1],zbin * DetectorReso[2]);
                ThreeVector<float> Dxyz = {(float) Dx->GetBinContent(xbin+1,ybin+1,zbin+1),(float) Dy->GetBinContent(xbin+1,ybin+1,zbin+1),(float) Dz->GetBinContent(xbin+1,ybin+1,zbin+1)};
                ThreeVector<float> True = RecoGrid + Dxyz;

                ThreeVector<float> RecoGrid_next((xbin+1)*DetectorReso[0],ybin*DetectorReso[1]+TPCVolume.GetDetectorOffset()[1],zbin*DetectorReso[2]);
                ThreeVector<float> Dxyz_next = {(float) Dx->GetBinContent(xbin+2,ybin+1,zbin+1),(float) Dy->GetBinContent(xbin+2,ybin+1,zbin+1),(float) Dz->GetBinContent(xbin+2,ybin+1,zbin+1)};
                ThreeVector<float> True_next = RecoGrid_next + Dxyz_next;

                ThreeVector<float> Rn = True_next - True;
                Position[xbin+ybin*(Resolution[0]-1)+zbin*(Resolution[0]-1)*Resolution[1]] = True + (float) 0.5 * Rn;
            }
        }
    }
    return Position;
}

// Write Emap into TH3 and store in root file
void WriteEmapRoot(std::vector<ThreeVector<float>>& Efield, TPCVolumeHandler& TPCVolume)
{
    // Store TPC properties which are important for the TH3 generation
    ThreeVector<unsigned long> Resolution = TPCVolume.GetDetectorResolution();
    ThreeVector<float> MinimumCoord = TPCVolume.GetMapMinimum();
    ThreeVector<float> MaximumCoord = TPCVolume.GetMapMaximum();

    // Initialize all TH3F
    std::vector<TH3F> Emap;
    Emap.push_back(TH3F("Emap_X","E field map X",Resolution[0],MinimumCoord[0],MaximumCoord[0],Resolution[1],MinimumCoord[1],MaximumCoord[1],Resolution[2],MinimumCoord[2],MaximumCoord[2]));
    Emap.push_back(TH3F("Emap_Y","E field map Y",Resolution[0],MinimumCoord[0],MaximumCoord[0],Resolution[1],MinimumCoord[1],MaximumCoord[1],Resolution[2],MinimumCoord[2],MaximumCoord[2]));
    Emap.push_back(TH3F("Emap_Z","E field map Z",Resolution[0],MinimumCoord[0],MaximumCoord[0],Resolution[1],MinimumCoord[1],MaximumCoord[1],Resolution[2],MinimumCoord[2],MaximumCoord[2]));

    // the loop should be consistent to the one in the EInterpolateMap()
    for(unsigned xbin = 0; xbin < Resolution[0]; xbin++)
    {
        for(unsigned ybin = 0; ybin < Resolution[1]; ybin++)
        {
            for(unsigned zbin = 0; zbin < Resolution[2]; zbin++)
            {
                // Loop over all coordinates dx,dy,dz
                for(unsigned coord = 0; coord < 3; coord++)
                {
                    // Fill interpolated grid points into histograms. bin=0 is underflow, bin = nbin+1 is overflow
                    Emap[coord].SetBinContent(xbin+1,ybin+1,zbin+1, Efield[zbin+ybin*Resolution[2]+xbin*Resolution[2]*Resolution[1]][coord]);
                } // end coordinate loop
            } // end zbin loop
        } // end ybin loop
    } // end zbin loop

    // Open and recreate output file
    TFile OutputFile("Emap.root", "recreate");

    // Loop over space coordinates
    for(unsigned coord = 0; coord < Emap.size(); coord++)
    {
        // Write every TH3 map into file
        Emap[coord].Write();
    }

    // Close output file and clean up
    OutputFile.Close();
    gDirectory->GetList()->Delete();
}