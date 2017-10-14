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
// *root_name is the name of the correction map
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

// This function uses Faraday's Law to calculate Ex at the boundary where the mesh didn't cover
// This function uses Maxwell-Faraday equation (Faraday's law of induction) that Ey, Ez = 0 at anode and cathode
// In principle, the methode can be generalised to extend E field in more generalised location of TPC,
// but we will lose the condition at the induction...
// 0 + -V_Lx + Int[(0 + Ez_dist(x,y,z))*dz,0,Lz] + Int[(E_drift + E_dist(x,y,z))*dx,0,Lx] = 0
// Int[(Ez_dist(x,y,z))*dz,0,Lz] + Int[(E_dist(x,y,z))*dx,0,Lx] = 0
// Coord is (xbin, ybin, zbin)
float EdgeEx(std::vector<ThreeVector<float>> &Efield, ThreeVector<unsigned long> Resolution, ThreeVector<float> Coord){


    for (unsigned xbin = 0; xbin < Resolution[0]; xbin++) {
        for (unsigned ybin = 0; ybin < Resolution[1]; ybin++) {
            for (unsigned zbin = 0; zbin < Resolution[2]; zbin++) {
                // Loop over all coordinates dx,dy,dz
                for (unsigned coord = 0; coord < 3; coord++) {
                    // Fill interpolated grid points into histograms. bin=0 is underflow, bin = nbin+1 is overflow
                    Emap[coord].SetBinContent(xbin + 1, ybin + 1, zbin + 1, Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][coord]);
                } // end coordinate loop
//                std::cout<<"xbin: "<<xbin<<"; ybin: "<<ybin<<"; zbin: "<<zbin<<"---Ex: "<<Efield[zbin+ybin*Resolution[2]+xbin*Resolution[2]*Resolution[1]][0]<<"; Ey: "<<Efield[zbin+ybin*Resolution[2]+xbin*Resolution[2]*Resolution[1]][1]<<"; Ez: "<< Efield[zbin+ybin*Resolution[2]+xbin*Resolution[2]*Resolution[1]][2]<<std::endl;
            } // end zbin loop
        } // end ybin loop
    } // end zbin loop

}



