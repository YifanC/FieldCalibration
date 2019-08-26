//
// Created by Yifan Chen on 08.10.17.
//
#include <vector>
#include <iostream>
#include <numeric>
#include "TROOT.h"
#include "TFile.h"
#include "TH3.h"
#include "../include/ThreeVector.hpp"
#include "../include/TPCVolumeHandler.hpp"
#include "../include/DriftVelocity.hpp"


// pair<En,Position>
// *root_name is the name of the correction map (input of the calculation)
std::pair<std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>>
Efield(TPCVolumeHandler &TPCVolume, float cryoTemp, float E0, float v0, const char *root_name) {

    float float_max = std::numeric_limits<float>::max();
    ThreeVector<float> Unknown = {float_max, float_max, float_max};

    TFile *InFile = new TFile(root_name, "READ");

//    TH3F *Dx = (TH3F *) InFile->Get("Reco_Displacement_X");
//    TH3F *Dy = (TH3F *) InFile->Get("Reco_Displacement_Y");
//    TH3F *Dz = (TH3F *) InFile->Get("Reco_Displacement_Z");

    TH3F *Dx = (TH3F *) InFile->Get("RecoBkwd_Displacement_X_Pos");
    TH3F *Dy = (TH3F *) InFile->Get("RecoBkwd_Displacement_Y_Pos");
    TH3F *Dz = (TH3F *) InFile->Get("RecoBkwd_Displacement_Z_Pos");


    ThreeVector<unsigned long> NrGrid = TPCVolume.GetDetectorResolution();
    ThreeVector<float> DetectorReso = {TPCVolume.GetDetectorSize()[0] / (NrGrid[0] - 1),
                                       TPCVolume.GetDetectorSize()[1] / (NrGrid[1] - 1),
                                       TPCVolume.GetDetectorSize()[2] / (NrGrid[2] - 1)};
    float Delta_x = DetectorReso[0]; //cm

    std::vector<ThreeVector<float>> En;
    std::vector<ThreeVector<float>> Position;

    for (unsigned Nz = 0; Nz < NrGrid[2]; Nz++) {
        for (unsigned Ny = 0; Ny < NrGrid[1]; Ny++) {
            // Since the E field calculation is based on the gap of the D map,
            // the number of x loop is one less than Resolution of the displacement map
//            for (unsigned Nx = 0; Nx < (NrGrid[0] - 1); Nx++) {
            for (unsigned Nx = NrGrid[0]; Nx > 1; Nx--) {
                //x = 0 (anode); x = Nx (cathode)
                ThreeVector<float> RecoGrid(Nx * DetectorReso[0],
                                            Ny * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
                                            Nz * DetectorReso[2]);
                ThreeVector<float> Dxyz = {(float) Dx->GetBinContent(Nx + 1, Ny + 1, Nz + 1),
                                           (float) Dy->GetBinContent(Nx + 1, Ny + 1, Nz + 1),
                                           (float) Dz->GetBinContent(Nx + 1, Ny + 1, Nz + 1)};

                if(Dxyz == Unknown) continue;

                ThreeVector<float> True = RecoGrid + Dxyz;

                ThreeVector<float> RecoGrid_next((Nx - 1) * DetectorReso[0],
                                                 Ny * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
                                                 Nz * DetectorReso[2]);
//                ThreeVector<float> RecoGrid_next((Nx + 1) * DetectorReso[0],
//                                                 Ny * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
//                                                 Nz * DetectorReso[2]);
                ThreeVector<float> Dxyz_next = {(float) Dx->GetBinContent(Nx, Ny + 1, Nz + 1),
                                                (float) Dy->GetBinContent(Nx, Ny + 1, Nz + 1),
                                                (float) Dz->GetBinContent(Nx, Ny + 1, Nz + 1)};
//                ThreeVector<float> Dxyz_next = {(float) Dx->GetBinContent(Nx + 2, Ny + 1, Nz + 1),
//                                                (float) Dy->GetBinContent(Nx + 2, Ny + 1, Nz + 1),
//                                                (float) Dz->GetBinContent(Nx + 2, Ny + 1, Nz + 1)};

                if(Dxyz_next == Unknown) continue;

                ThreeVector<float> True_next = RecoGrid_next + Dxyz_next;

                ThreeVector<float> Rn = True_next - True;

                // mm/us, the magnitude of the drift velocity at the local gap
                // Be very careful that Rn.GetNorm() and Delta_x are in the unit of cm, not mm
                float vn = Rn.GetNorm() / Delta_x * v0;

                // Only include the valid calculated Efield (which is not float_max) in the pre-mesh
                if (searchE(vn, cryoTemp, E0) < 0.5 * std::numeric_limits<float>::max()) {
                    // the E field as a vector has the same direction of Rn (vector) in each gap
                    // E kV/cm
                    En.push_back(searchE(vn, cryoTemp, E0) / Rn.GetNorm() * Rn);
                    // Set the local E field (E field from each gap) at the middle of the gap
                    Position.push_back(True + (float) 0.5 * Rn);
                }
            }
        }
    }

    std::pair<std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>> field;
    field = std::make_pair(En,Position);

    return field;
}

// Same as Efield, just take vector as input instead of the root file
// tuple<vn, En, Position>
std::tuple<std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>>
EfieldvecMap(TPCVolumeHandler &TPCVolume, float cryoTemp, float E0, float v0, std::vector<ThreeVector<float>> DMapTT) {

//    TFile *InFile = new TFile(root_name, "READ");
//
//    TH3F *Dx = (TH3F *) InFile->Get("Reco_Displacement_X");
//    TH3F *Dy = (TH3F *) InFile->Get("Reco_Displacement_Y");
//    TH3F *Dz = (TH3F *) InFile->Get("Reco_Displacement_Z");

    float float_max = std::numeric_limits<float>::max();
    ThreeVector<float> Unknown = {float_max, float_max, float_max};


    ThreeVector<unsigned long> NrGrid = TPCVolume.GetDetectorResolution();
    ThreeVector<float> DetectorReso = {TPCVolume.GetDetectorSize()[0] / (NrGrid[0] - 1),
                                       TPCVolume.GetDetectorSize()[1] / (NrGrid[1] - 1),
                                       TPCVolume.GetDetectorSize()[2] / (NrGrid[2] - 1)};
    float Delta_x = DetectorReso[0]; //cm

    std::vector<ThreeVector<float>> Velocity;
    std::vector<ThreeVector<float>> En;
    std::vector<ThreeVector<float>> Position;

    for (unsigned Nz = 0; Nz < NrGrid[2]; Nz++) {
        for (unsigned Ny = 0; Ny < NrGrid[1]; Ny++) {
            // To save the En, vn from last calculation step, initiate the containers and the flags
            ThreeVector<float> Elast;
            bool AvgLast = false;
            bool AvgThis = false;

            // Since the E field calculation is based on the gap of the D map,
            // the number of x loop is one less than Resolution of the displacement map
//            for (unsigned Nx = 0; Nx < (NrGrid[0] - 1); Nx++) {
            for (unsigned Nx = NrGrid[0]; Nx > 1; Nx--) {
                //x = 0 (anode); x = Nx (cathode)
                ThreeVector<float> RecoGrid(Nx * DetectorReso[0] + TPCVolume.GetDetectorOffset()[0],
                                            Ny * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
                                            Nz * DetectorReso[2] + TPCVolume.GetDetectorOffset()[2]);

//                ThreeVector<float> Dxyz = {(float) Dx->GetBinContent(Nx + 1, Ny + 1, Nz + 1),
//                                           (float) Dy->GetBinContent(Nx + 1, Ny + 1, Nz + 1),
//                                           (float) Dz->GetBinContent(Nx + 1, Ny + 1, Nz + 1)};
                ThreeVector<float> Dxyz = DMapTT[Nz + (NrGrid[2] * (Ny + NrGrid[1] * Nx))];

                if(Dxyz == Unknown) continue;

                ThreeVector<float> True = RecoGrid + Dxyz;

//                ThreeVector<float> RecoGrid_next((Nx + 1) * DetectorReso[0] + TPCVolume.GetDetectorOffset()[0],
//                                                 Ny * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
//                                                 Nz * DetectorReso[2] + TPCVolume.GetDetectorOffset()[2]);
                ThreeVector<float> RecoGrid_next((Nx - 1) * DetectorReso[0] + TPCVolume.GetDetectorOffset()[0],
                                                 Ny * DetectorReso[1] + TPCVolume.GetDetectorOffset()[1],
                                                 Nz * DetectorReso[2] + TPCVolume.GetDetectorOffset()[2]);

//                ThreeVector<float> Dxyz_next = {(float) Dx->GetBinContent(Nx + 2, Ny + 1, Nz + 1),
//                                                (float) Dy->GetBinContent(Nx + 2, Ny + 1, Nz + 1),
//                                                (float) Dz->GetBinContent(Nx + 2, Ny + 1, Nz + 1)};

//                ThreeVector<float> Dxyz_next = DMapTT[Nz + (NrGrid[2] * (Ny + NrGrid[1] * (Nx+1)))];
                ThreeVector<float> Dxyz_next = DMapTT[Nz + (NrGrid[2] * (Ny + NrGrid[1] * (Nx - 1)))];

                if(Dxyz_next == Unknown) continue;

                ThreeVector<float> True_next = RecoGrid_next + Dxyz_next;

                ThreeVector<float> Rn = True_next - True;

                // mm/us, the magnitude of the drift velocity at the local gap
                // Be very careful that Rn.GetNorm() and Delta_x are in the unit of cm, not mm
                float vn = Rn.GetNorm() / Delta_x * v0;
                AvgThis = false;

                // Only include the valid calculated Efield (which is not float_max) in the pre-mesh
                if (searchE(vn, cryoTemp, E0) < 0.5 * std::numeric_limits<float>::max()) {

                    AvgThis = true;

                    // the drift velocity as a vector has the same direction of Rn (vector) in each gap
                    // v mm/us
                    Velocity.push_back(vn / Rn.GetNorm() * Rn);
                    // the E field as a vector has the same direction of Rn (vector) in each gap
                    // E kV/cm
                    ThreeVector<float> Ethis = searchE(vn, cryoTemp, E0) / Rn.GetNorm() * Rn;
                    En.push_back(Ethis);
                    // Set the local E field (E field from each gap) at the middle of the gap
                    Position.push_back(True + (float) 0.5 * Rn);

                    std::cout<<"A"<<std::endl;

//                    if(Nx == 0){
//                        Velocity.push_back(vn / Rn.GetNorm() * Rn);
//                        En.push_back(searchE(vn, cryoTemp, E0) / Rn.GetNorm() * Rn);
//                        Position.push_back(True);
//                    }
//
//                    //Fill the boundary as inner value
//                    if(Nx == NrGrid[0] - 2){
//                        Velocity.push_back(vn / Rn.GetNorm() * Rn);
//                        En.push_back(searchE(vn, cryoTemp, E0) / Rn.GetNorm() * Rn);
//                        Position.push_back(True_next);
//                    }

                    if(Nx == NrGrid[0]){
                        Velocity.push_back(vn / Rn.GetNorm() * Rn);
                        En.push_back(searchE(vn, cryoTemp, E0) / Rn.GetNorm() * Rn);
                        Position.push_back(True);
                    }

                    //Fill the boundary as inner value
                    if(Nx == 2) {
                        Velocity.push_back(vn / Rn.GetNorm() * Rn);
                        En.push_back(searchE(vn, cryoTemp, E0) / Rn.GetNorm() * Rn);
                        Position.push_back(True_next);
                    }

                    //If the previous step has invalid value or Nx==0, just fill the joint with result from the current step
                    if(!AvgLast){
                        Velocity.push_back(vn / Rn.GetNorm() * Rn);
                        En.push_back(Ethis);
                        Position.push_back(True);
                    }
                    //Fill the joint with average E-field and corresponding drift v
                    if(Nx > 0 && AvgLast){
                        // If AvgLast is true, Elast exists with correct value
                        ThreeVector<float> Eaverage = (float)0.5 * (Ethis + Elast);
                        Velocity.push_back(ElectronDriftVelocity(cryoTemp, Eaverage.GetNorm()) / Eaverage.GetNorm() * Eaverage);
                        En.push_back(Eaverage);
                        Position.push_back(True);

                    }
                    std::cout<<Nx << "B"<<std::endl;
                    // store E-field vector for next calculation step
                    Elast = searchE(vn, cryoTemp, E0) / Rn.GetNorm() * Rn;
                }
                AvgLast = AvgThis;
            }
        }
    }

//    std::pair<std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>> field;
//    field = std::make_pair(En,Position);
    std::cout<<"C"<<std::endl;
    auto vnE = std::make_tuple(Velocity, En, Position);

    return vnE;
}


// tuple <Ex, Ey, Ez, PositionX, PositionY, PositionZ>
// *root_name is the name of the correction map (input of the calculation)
// The boundary condition at TPC edge is included in this function by default.
std::tuple<std::vector<float >, std::vector<float >, std::vector<float >,
        std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>>
EfieldXYZwithBoundary(TPCVolumeHandler &TPCVolume, float cryoTemp, float E0, float v0, const char *root_name) {

    TFile *InFile = new TFile(root_name, "READ");

    TH3F *Dx = (TH3F *) InFile->Get("Reco_Displacement_X");
    TH3F *Dy = (TH3F *) InFile->Get("Reco_Displacement_Y");
    TH3F *Dz = (TH3F *) InFile->Get("Reco_Displacement_Z");

    ThreeVector<unsigned long> NrGrid = TPCVolume.GetDetectorResolution();
    ThreeVector<float> DetectorReso = {TPCVolume.GetDetectorSize()[0] / (NrGrid[0] - 1),
                                       TPCVolume.GetDetectorSize()[1] / (NrGrid[1] - 1),
                                       TPCVolume.GetDetectorSize()[2] / (NrGrid[2] - 1)};
    float Delta_x = DetectorReso[0]; //cm
    ThreeVector<float> Offset = TPCVolume.GetDetectorOffset();

    ThreeVector<float> MinimumCoord = TPCVolume.GetMapMinimum();
    ThreeVector<float> MaximumCoord = TPCVolume.GetMapMaximum();

    std::vector<float> Ex;
    std::vector<float> Ey;
    std::vector<float> Ez;
    std::vector<ThreeVector<float>> Position;
    std::vector<ThreeVector<float>> PositionX;
    std::vector<ThreeVector<float>> PositionY;
    std::vector<ThreeVector<float>> PositionZ;

    for (unsigned Nz = 0; Nz < NrGrid[2]; Nz++) {
        for (unsigned Ny = 0; Ny < NrGrid[1]; Ny++) {
            // Since the E field calculation is based on the gap of the D map,
            // the number of x loop is one less than Resolution of the displacement map
            for (unsigned Nx = 0; Nx < (NrGrid[0] - 1); Nx++) {
                //x = 0 (anode); x = NrGrid[0] (cathode)
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
                    // E kV/cm
//                    En.push_back(searchE(vn, cryoTemp, E0) / Rn.GetNorm() * Rn);
                    ThreeVector<float> En = searchE(vn, cryoTemp, E0) / Rn.GetNorm() * Rn;
                    Ex.push_back(En[0]);
                    Ey.push_back(En[1]);
                    Ez.push_back(En[2]);
                    // Set the local E field (E field from each gap) at the middle of the gap
                    Position.push_back(True + (float) 0.5 * Rn);
                }


            }
        }
    }

    PositionX = Position;
    PositionY = Position;
    PositionZ = Position;

    // The following part can be merge into the former big loop.
    // However, it may be more clear in this way.
    // If you would like to merge into the previous loop, be cautious the mesh at Xbin = NrGrid[0]-1.
    // It's better to have larger mesh size.
    // X=Xmin, X=Xmax
    for (unsigned ybin = 0; ybin < NrGrid[1]; ybin++) {
        for (unsigned zbin = 0; zbin < NrGrid[2]; zbin++) {

            //Push back the location (x,y,z coord) of Anode Points. Anode sits at x=0
            ThreeVector<float> gridX0 = {MinimumCoord[0], ybin * DetectorReso[1] + Offset[1], zbin * DetectorReso[2] + Offset[2]};
            PositionX.push_back(gridX0);
            Ey.push_back(0.);
            Ez.push_back(0.);

            //Push back the location (x,y,z coord) of Anode Points. Anode sits at x=0
            ThreeVector<float> gridX1 = {MaximumCoord[0], ybin * DetectorReso[1] + Offset[1], zbin * DetectorReso[2] + Offset[2]};
            PositionX.push_back(gridX1);
            Ey.push_back(0.);
            Ez.push_back(0.);

        }
    }

//    // Y=Ymin, Y=Ymax
//    for (unsigned xbin = 0; xbin < NrGrid[0]; xbin++) {
//        for (unsigned zbin = 0; zbin < NrGrid[2]; zbin++) {
//
//            //Push back the location (x,y,z coord) of Anode Points. Anode sits at x=0
//            ThreeVector<float> gridY0 = {xbin * DetectorReso[0] + Offset[0], MinimumCoord[1], zbin * DetectorReso[2] + Offset[2]};
//            PositionY.push_back(gridY0);
//            Ex.push_back(0.);
//            Ez.push_back(0.);
//
//            //Push back the location (x,y,z coord) of Anode Points. Anode sits at x=0
//            ThreeVector<float> gridY1 = {xbin * DetectorReso[0] + Offset[0], MaximumCoord[1], zbin * DetectorReso[2] + Offset[2]};
//            PositionY.push_back(gridY1);
//            Ex.push_back(0.);
//            Ez.push_back(0.);
//
//        }
//    }
//
//    // Z=Zmin, Z=Zmax
//    for (unsigned ybin = 0; ybin < NrGrid[1]; ybin++) {
//        for (unsigned xbin = 0; xbin < NrGrid[0]; xbin++) {
//
//            //Push back the location (x,y,z coord) of Anode Points. Anode sits at x=0
//            ThreeVector<float> gridZ0 = {xbin * DetectorReso[0] + Offset[0], ybin * DetectorReso[1] + Offset[1], MinimumCoord[2]};
//            PositionZ.push_back(gridZ0);
//            Ey.push_back(0.);
//            Ex.push_back(0.);
//
//            //Push back the location (x,y,z coord) of Anode Points. Anode sits at x=0
//            ThreeVector<float> gridZM = {xbin * DetectorReso[0] + Offset[0], ybin * DetectorReso[1] + Offset[1], MaximumCoord[2]};
//            PositionZ.push_back(gridZM);
//            Ey.push_back(0.);
//            Ex.push_back(0.);
//
//        }
//    }

    std::cout<<"Position X size: "<<PositionX.size()
            <<"Position Y size: "<<PositionY.size()
            <<"Position Z size: "<<PositionZ.size()
                                 <<std::endl;
    std::tuple<std::vector<float >, std::vector<float >, std::vector<float >,
            std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>> field;
    field = std::make_tuple(Ex, Ey, Ez, PositionX, PositionY, PositionZ);

    return field;
}



// This function uses Faraday's Law to calculate Ex at the boundary where the mesh didn't cover
// This function uses Maxwell-Faraday equation (Faraday's law of induction) that Ey, Ez = 0 at anode and cathode
// In principle, the methode can be generalised to extend E field in more generalised location of TPC,
// but we will lose the condition at the induction...
// We choose Integral path along z instead of y because of smaller uncertainty
// 0 + -V_Lx + Int[(0 + Ez_dist(x,y,z))*dz,0,Lz] + Int[(E_drift + E_dist(x,y,z))*dx,0,Lx] = 0
// Int[(Ez_dist(x,y,z))*dz,0,Lz] + Int[(E_dist(x,y,z))*dx,0,Lx] = 0
// Coord is (xbin, ybin, zbin) from 0 to Resolution-1
// return the distortion in Ex
float EdgeEx(std::vector<ThreeVector<float>> &Efield, ThreeVector<unsigned long> Resolution, ThreeVector<float> Unit, float E0, ThreeVector<unsigned long> Coord){
    std::vector<float> Ex;
    std::vector<float> vec_Int_yz;
    float Integral_yz;
    float Integral_x=0;


    if(Coord[0]!=0 && Coord[0]!=Resolution[0]-1){
        std::cerr<<"We only calculate Ex distortion on Anode or Cathode now"<<std::endl;
    }


//    for(int nx = 0+1;nx<Resolution[0]-1;nx++){
    for(int i = -1;i<2;i++){
        unsigned long Xend = Resolution[0]*(0.5+i*0.3);
//        unsigned long Xend = Resolution[0]-2;
//        unsigned long Xend = nx;
        Integral_x =0;
//        std::cout<<"Xend: "<<Xend<<std::endl;

        //Integral path always starts from (x,y,z) to (Xend, y,z) as first leg
        //if anode
        if(Coord[0]==0){
            for (int xbin = 0+1; xbin < Xend+1; xbin++) {
                if(xbin == Xend){
                    Integral_x += (Efield[Coord[2] + Coord[1] * Resolution[2] + xbin * Resolution[2] * Resolution[1]][0] - E0) * (Unit[0]*0.5);
                }
                else{
                    Integral_x += (Efield[Coord[2] + Coord[1] * Resolution[2] + xbin * Resolution[2] * Resolution[1]][0] - E0) * Unit[0];
                }
//                std::cout<<"xbin: "<<xbin<<"... sum Integral x: "<<Integral_x<<std::endl;
            }
        }

        //if cathode
        if(Coord[0]==Resolution[0]-1){
            for (int xbin = Resolution[0]-1; xbin > Xend-1; xbin--) {
                if(xbin == Xend){
                    Integral_x -= (Efield[Coord[2] + Coord[1] * Resolution[2] + xbin * Resolution[2] * Resolution[1]][0] - E0) * (Unit[0]*0.5);
                }
                else{
                    Integral_x -= (Efield[Coord[2] + Coord[1] * Resolution[2] + xbin * Resolution[2] * Resolution[1]][0] - E0) * Unit[0];
                }
            }
        }

//        std::cout<<"X: Integral_x: "<<Integral_x<<std::endl;


////////////////////
        Integral_yz =0;
        for (int zbin = Coord[2]; zbin < Resolution[2]; zbin++) {
//            std::cout<<"zbin: "<<zbin<<std::endl;
            //The direction is + here
            if((zbin == Coord[2]) || (zbin == Resolution[2] - 1)){
                Integral_yz += Efield[zbin + Coord[1] * Resolution[2] + Xend * Resolution[2] * Resolution[1]][2]*(Unit[2]*0.5);
            }
            else{
                Integral_yz += Efield[zbin + Coord[1] * Resolution[2] + Xend * Resolution[2] * Resolution[1]][2]*Unit[2];
            }

        }
//        std::cout<<"Integral_yz: "<<Integral_yz<<std::endl;

        if(Coord[0]==0){Ex.push_back((0 - Integral_x - Integral_yz)/ (Unit[0]*0.5));}
        if(Coord[0]==Resolution[0]-1){Ex.push_back((Integral_x + Integral_yz)/ (Unit[0]*0.5));}

//        std::cout<<"xbin: "<<nx<<", Integral_x: "<<Integral_x<<", Integral_yz: "<<Integral_yz<<", Ex: "<<Ex.back()<<", percentage : "<<Ex.back()/E0*100<<std::endl;
///////////////////////////
        Integral_yz =0;
        for (int zbin = Coord[2]; zbin > -1; zbin--) {
            //The direction is + here
            if((zbin == Coord[2]) || (zbin == 0)){
                Integral_yz -= Efield[zbin + Coord[1] * Resolution[2] + Xend * Resolution[2] * Resolution[1]][2]*(Unit[2]*0.5);
            }
            else{
                Integral_yz -= Efield[zbin + Coord[1] * Resolution[2] + Xend * Resolution[2] * Resolution[1]][2]*Unit[2];

            }

        }

        if(Coord[0]==0){Ex.push_back((0 - Integral_x - Integral_yz)/ (Unit[0]*0.5));}
        if(Coord[0]==Resolution[0]-1){Ex.push_back((Integral_x + Integral_yz)/ (Unit[0]*0.5));}

        /////////////////////////
        Integral_yz =0;
        for (int ybin = Coord[1]; ybin < Resolution[1]; ybin++) {

            //The direction is + here
            if((ybin == Coord[1]) || (ybin == Resolution[1]-1)){
                Integral_yz += Efield[Coord[2] + ybin * Resolution[2] + Xend * Resolution[2] * Resolution[1]][1]*(Unit[1]*0.5);
            }
            else{
                Integral_yz += Efield[Coord[2] + ybin * Resolution[2] + Xend * Resolution[2] * Resolution[1]][1]*Unit[1];
            }

        }

        if(Coord[0]==0){Ex.push_back((0 - Integral_x - Integral_yz)/ (Unit[0]*0.5));}
        if(Coord[0]==Resolution[0]-1){Ex.push_back((Integral_x + Integral_yz)/ (Unit[0]*0.5));}

        /////////////////////////
        Integral_yz =0;
        for (int ybin = Coord[1]; ybin > -1; ybin--) {
            //The direction is + here
            if((ybin == Coord[1]) || (ybin == 0)){
                Integral_yz -= Efield[Coord[2] + ybin * Resolution[2] + Xend * Resolution[2] * Resolution[1]][1]*(Unit[1]*0.5);
            }
            else{
                Integral_yz -= Efield[Coord[2] + ybin * Resolution[2] + Xend * Resolution[2] * Resolution[1]][1]*Unit[1];
            }

        }

        if(Coord[0]==0){Ex.push_back((0 - Integral_x - Integral_yz)/ (Unit[0]*0.5));}
        if(Coord[0]==Resolution[0]-1){Ex.push_back((Integral_x + Integral_yz)/ (Unit[0]*0.5));}
    }


    float meanDistEx = std::accumulate(Ex.begin(), Ex.end(),0.0)/Ex.size();

    return meanDistEx+E0;
}


