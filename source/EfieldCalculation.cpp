//
// Created by Yifan Chen on 08.10.17.
//
#include "TROOT.h"
#include "TFile.h"
#include "TH3.h"
#include "../include/ThreeVector.hpp"
#include "../include/TPCVolumeHandler.hpp"
#include "../include/DriftVelocity.hpp"

// pair<En,Position>
std::pair<std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>>
Efield(TPCVolumeHandler &TPCVolume, float cryoTemp, float E0, float v0, const char *root_name) {

    TFile *InFile = new TFile(root_name, "READ");

    TH3F *Dx = (TH3F *) InFile->Get("Reco_Displacement_X");
    TH3F *Dy = (TH3F *) InFile->Get("Reco_Displacement_Y");
    TH3F *Dz = (TH3F *) InFile->Get("Reco_Displacement_Z");

    ThreeVector<unsigned long> NrGrid = TPCVolume.GetDetectorResolution();
    ThreeVector<float> DetectorReso = {TPCVolume.GetDetectorSize()[0] / (NrGrid[0] - 1),
                                       TPCVolume.GetDetectorSize()[1] / (NrGrid[1] - 1),
                                       TPCVolume.GetDetectorSize()[2] / (NrGrid[2] - 1)};
    float Delta_x = DetectorReso[0]; //cm

//    std::vector<ThreeVector< float>> En(Resolution[2] * Resolution[1] * (Resolution[0]-1));
//    std::vector<ThreeVector<float>> Position(Resolution[2]*Resolution[1]*(Resolution[0]-1));
    std::vector<ThreeVector<float>> En;
    std::vector<ThreeVector<float>> Position;

    for (unsigned Nz = 0; Nz < NrGrid[2]; Nz++) {
        for (unsigned Ny = 0; Ny < NrGrid[1]; Ny++) {
            // Since the E field calculation is based on the gap of the D map, the number of x loop is one less than Resolution of the displacement map
            for (unsigned Nx = 0; Nx < (NrGrid[0] - 1); Nx++) {
                //x = 0 (anode); x = Nx (cathode)
                ThreeVector<float> RecoGrid(Nx * DetectorReso[0],
                                            Ny * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
                                            Nz * DetectorReso[2]);
                ThreeVector<float> Dxyz = {(float) Dx->GetBinContent(Nx + 1, Ny + 1, Nz + 1),
                                           (float) Dy->GetBinContent(Nx + 1, Ny + 1, Nz + 1),
                                           (float) Dz->GetBinContent(Nx + 1, Ny + 1, Nz + 1)};
                ThreeVector<float> True = RecoGrid + Dxyz;

                ThreeVector<float> RecoGrid_next((Nx + 1) * DetectorReso[0],
                                                 Ny * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
                                                 Nz * DetectorReso[2]);
                ThreeVector<float> Dxyz_next = {(float) Dx->GetBinContent(Nx + 2, Ny + 1, Nz + 1),
                                                (float) Dy->GetBinContent(Nx + 2, Ny + 1, Nz + 1),
                                                (float) Dz->GetBinContent(Nx + 2, Ny + 1, Nz + 1)};
                ThreeVector<float> True_next = RecoGrid_next + Dxyz_next;

                ThreeVector<float> Rn = True_next - True;

                // mm/us, the magnitude of the drift velocity at the local gap
                // Be very careful that Rn.GetNorm() and Delta_x are in the unit of cm, not mm
                float vn = Rn.GetNorm() / Delta_x * v0;


                if (searchE(vn, cryoTemp, E0) < 0.5 * std::numeric_limits<float>::max()) {
                    // the E field as a vector has the same direction of Rn (vector) in each gap
                    En.push_back(searchE(vn, cryoTemp, E0) / Rn.GetNorm() * Rn);
                    // Set the local E field (E field from each gap) at the middle of the gap
                    Position.push_back(True + (float) 0.5 * Rn);
                }
//                Position[Nx+ybin*(Resolution[0]-1)+zbin*(Resolution[0]-1)*Resolution[1]] = True + (float) 0.5 * Rn;
//                En[Nx+ybin*(Resolution[0]-1)+zbin*(Resolution[0]-1)*Resolution[1]] = searchE(vn,cryoTemp,E0) / Rn.GetNorm() * Rn;
            }
        }
    }

    std::pair<std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>> field;
    field = std::make_pair(En,Position);

    return field;
}


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
//            for (unsigned Nx = 0; Nx < (Resolution[0] - 1); Nx++) {
//                //Nx =1, x=5 close to the anode; Nx = Nx, x=250 close to the cathode
//                ThreeVector<float> RecoGrid(Nx * DetectorReso[0],
//                                            ybin * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
//                                            zbin * DetectorReso[2]);
//                ThreeVector<float> Dxyz = {(float) Dx->GetBinContent(Nx + 1, ybin + 1, zbin + 1),
//                                           (float) Dy->GetBinContent(Nx + 1, ybin + 1, zbin + 1),
//                                           (float) Dz->GetBinContent(Nx + 1, ybin + 1, zbin + 1)};
//                ThreeVector<float> True = RecoGrid + Dxyz;
//
//                ThreeVector<float> RecoGrid_next((Nx + 1) * DetectorReso[0],
//                                                 ybin * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
//                                                 zbin * DetectorReso[2]);
//                ThreeVector<float> Dxyz_next = {(float) Dx->GetBinContent(Nx + 2, ybin + 1, zbin + 1),
//                                                (float) Dy->GetBinContent(Nx + 2, ybin + 1, zbin + 1),
//                                                (float) Dz->GetBinContent(Nx + 2, ybin + 1, zbin + 1)};
//                ThreeVector<float> True_next = RecoGrid_next + Dxyz_next;
//
//                ThreeVector<float> Rn = True_next - True;
//
//                // To synchronize the vector order of En and Position through the output of searchE (if the output is floatmax, then abandon both elements in Position and En)
//                float vn = Rn.GetNorm() / Delta_x * v0;
//                if (searchE(vn, cryoTemp, E0) < 0.5 * std::numeric_limits<float>::max()) {
//                    Position.push_back(searchE(vn, cryoTemp, E0) / Rn.GetNorm() * Rn);
//                }
////                Position[Nx+ybin*(Resolution[0]-1)+zbin*(Resolution[0]-1)*Resolution[1]] = True + (float) 0.5 * Rn;
//            }
//        }
//    }
//    std::cout << "Position size: " << Position.size() << std::endl;
//    return Position;
//}
//
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
//            for (unsigned Nx = 0; Nx < (TPCVolume.GetDetectorResolution()[0] - 1); Nx++) {
//                //Nx =1, x=5 close to the anode; Nx = Nx, x=250 close to the cathode
//                // the corner has the weight of the bin. Do not move the weight to the geometry center!
//                ThreeVector<float> RecoGrid(Nx * DetectorReso[0],
//                                            ybin * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
//                                            zbin * DetectorReso[2]);
//                ThreeVector<float> Dxyz = {(float) Dx->GetBinContent(Nx + 1, ybin + 1, zbin + 1),
//                                           (float) Dy->GetBinContent(Nx + 1, ybin + 1, zbin + 1),
//                                           (float) Dz->GetBinContent(Nx + 1, ybin + 1, zbin + 1)};
//                ThreeVector<float> True = RecoGrid + Dxyz;
//
//                ThreeVector<float> RecoGrid_next((Nx + 1) * DetectorReso[0],
//                                                 ybin * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
//                                                 zbin * DetectorReso[2]);
//                ThreeVector<float> Dxyz_next = {(float) Dx->GetBinContent(Nx + 2, ybin + 1, zbin + 1),
//                                                (float) Dy->GetBinContent(Nx + 2, ybin + 1, zbin + 1),
//                                                (float) Dz->GetBinContent(Nx + 2, ybin + 1, zbin + 1)};
//                ThreeVector<float> True_next = RecoGrid_next + Dxyz_next;
//
//                ThreeVector<float> Rn = True_next - True;
//                float vn = Rn.GetNorm() / Delta_x * v0; // mm/us, the magnitude of the drift velocity at the local point
//                //To be very careful that Rn.GetNorm() and Delta_x are in the unit of cm, not mm
//
////                if(searchE(vn,cryoTemp,E0)>0.6){
////                    std::cout<<"Need investigation! Nx: "<<Nx<<"; ybin: "<<ybin<<"; zbin: "<<zbin<<"; |Rn|: "<<Rn.GetNorm()<<"; |E|: "<<searchE(vn,cryoTemp,E0)<<std::endl;
////                }
//
//                // the E field as a vector has the same direction of Rn (Threevector)
//                if (searchE(vn, cryoTemp, E0) < 0.5 * std::numeric_limits<float>::max()) {
//                    En.push_back(searchE(vn, cryoTemp, E0) / Rn.GetNorm() * Rn);
//                }
//
////                En[Nx+ybin*(Resolution[0]-1)+zbin*(Resolution[0]-1)*Resolution[1]] = searchE(vn,cryoTemp,E0) / Rn.GetNorm() * Rn;
//            }
//        }
//    }
//
//    std::cout << "En size: " << En.size() << std::endl;
//    return En;
//}


