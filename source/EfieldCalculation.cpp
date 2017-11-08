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
//                std::cout<<"Rn.GetNorm: "<<Rn.GetNorm()<<", Delta_x: "<<Delta_x<<", v0: "<<v0<<std::endl;
//                std::cout<<"Rn_x: "<<Rn[0]<<", Rn_y: "<<Rn[1]<<", Rn_z: "<<Rn[2]<<std::endl;
//                std::cout<<"True_next: "<<True_next.GetNorm()<<", True: "<<True.GetNorm()<<std::endl;

                if (searchE(vn, cryoTemp, E0) < 0.5 * std::numeric_limits<float>::max()) {
                    // the E field as a vector has the same direction of Rn (vector) in each gap
                    // E kV/cm
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
//    unsigned long nx = Coord[0];
//    unsigned long ny = Coord[1];
//    unsigned long nz = Coord[2];
//    std::cout<<"Coord[0]: "<<Coord[0]<<", Coord[1]: "<<Coord[1]<<", Coord[2]: "<<Coord[2]<<std::endl;

    // In one xbin (xstream), at least one end should have max one Ex (on anode or cathode) missing.
//    if(Coord[0]!=0 && Coord[0]!=Resolution[0]-1){
//        std::cerr<<"We only calculate Ex distortion on Anode or Cathode now"<<std::endl;
//    }

    if(Coord[0]!=0 && Coord[0]!=Resolution[0]-1){
        std::cerr<<"We only calculate Ex distortion on Anode or Cathode now"<<std::endl;
    }

    //?????????
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
//                std::cout<<"Ey: "<< Efield[zbin + Coord[1] * Resolution[2] + Xend * Resolution[2] * Resolution[1]][2]<<std::endl;
            }

        }
//        std::cout<<"Integral_yz: "<<Integral_yz<<std::endl;

        if(Coord[0]==0){Ex.push_back((0 - Integral_x - Integral_yz)/ (Unit[0]*0.5));}
        if(Coord[0]==Resolution[0]-1){Ex.push_back((Integral_x + Integral_yz)/ (Unit[0]*0.5));}

        /////////////////////////
        Integral_yz =0;
        for (int ybin = Coord[1]; ybin < Resolution[1]; ybin++) {
//            std::cout<<"ybin: "<<ybin<<std::endl;
            //The direction is + here
            if((ybin == Coord[1]) || (ybin == Resolution[1]-1)){
                Integral_yz += Efield[Coord[2] + ybin * Resolution[2] + Xend * Resolution[2] * Resolution[1]][1]*(Unit[1]*0.5);
            }
            else{
                Integral_yz += Efield[Coord[2] + ybin * Resolution[2] + Xend * Resolution[2] * Resolution[1]][1]*Unit[1];
            }

        }
//        std::cout<<"Integral_yz: "<<Integral_yz<<std::endl;

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
//        std::cout<<"Integral_yz: "<<Integral_yz<<std::endl;

        if(Coord[0]==0){Ex.push_back((0 - Integral_x - Integral_yz)/ (Unit[0]*0.5));}
        if(Coord[0]==Resolution[0]-1){Ex.push_back((Integral_x + Integral_yz)/ (Unit[0]*0.5));}
    }

//    for(int n=0;n<Ex.size();n++){
//        std::cout<<"Ex: "<<Ex[n]<<", ";
//    }
//    std::cout<<std::endl;



    float meanDistEx = std::accumulate(Ex.begin(), Ex.end(),0.0)/Ex.size();
//    std::cout<<"mean Ditorted Ex: "<<meanDistEx+E0<<", (Ex-E0)/E0*100% : "<< meanDistEx/E0*100<< std::endl;

    return meanDistEx+E0;
}


/*
///////////////////////////////////////
std::vector<ThreeVector<float>>
MaxwellEmap(std::vector<ThreeVector<float>> &Efield, TPCVolumeHandler &TPCVolume,
            ThreeVector<unsigned long> Resolution, float E0){

    float float_max = std::numeric_limits<float>::max();

    // the loop should be consistent to the one in the EInterpolateMap()
//    for (unsigned xbin = 0; xbin < Resolution[0]; xbin++) {
    for (unsigned xbin = 0; xbin < 3; xbin++) {
        for (unsigned ybin = 0; ybin < Resolution[1]; ybin++) {
            for (unsigned zbin = 0; zbin < Resolution[2]; zbin++) {
//                if((abs(Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][0]) > 0.5*float_max ||
//                    abs(Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][1]) > 0.5*float_max ||
//                    abs(Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][2]) > 0.5*float_max))
//                {
                if(abs(Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][0]) > 0.5*float_max )
                {
                    if(abs(Efield[zbin + ybin * Resolution[2] + (xbin+1) * Resolution[2] * Resolution[1]][0]) < 0.5*float_max||
                       abs(Efield[zbin + ybin * Resolution[2] + (xbin-1) * Resolution[2] * Resolution[1]][0]) < 0.5*float_max)
                    {
                        //hope
                        //always two unknown vairables. take 2 path to solve 2 varibles
                    }
                    else
                    {
                        //what to do?????
                    }
                    Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][2];
                    ThreeVector<unsigned long> Coord = {xbin, ybin, zbin};
                    Emap[0].SetBinContent(xbin + 1, ybin + 1, zbin + 1, EdgeEx(Efield, Resolution, Unit, E0, Coord));
                    Emap[1].SetBinContent(xbin + 1, ybin + 1, zbin + 1, 0);
                    Emap[2].SetBinContent(xbin + 1, ybin + 1, zbin + 1, 0);
                    std::cout<<"xbin: "<<xbin<<", ybin: "<<ybin<<", zbin: "<<zbin<<", Ex distortion: "<<EdgeEx(Efield, Resolution, Unit, E0, Coord)<<std::endl;
                }
                if(abs(Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][1]) > 0.5*float_max ){

                    if(abs(Efield[zbin + (ybin+1) * Resolution[2] + xbin * Resolution[2] * Resolution[1]][1]) < 0.5*float_max||
                       abs(Efield[zbin + (ybin-1) * Resolution[2] + xbin * Resolution[2] * Resolution[1]][1]) < 0.5*float_max)
                    {
                        //hope
                        //maybe one unknown variable or 2.. safe to take two path
                        for (int xbin = X0; xbin < X1 + 1; xbin++){

                            if(xbin == X0 || xbin == X1){
                                Unit0 = DetectorReso[0]*0.5;
                            }
                            else{
                                Unit0 = DetectorReso[0];
                            }

                            Integral += (Ex->GetBinContent(xbin+1,Y0+1,Z0+1)-E0)*Unit0;
                            Integral -= (Ex->GetBinContent(xbin+1,Y1+1,Z1+1)-E0)*Unit0;

                            // Integral += Ex->GetBinContent(xbin+1,Y0+1,Z0+1)*Unit0;
                            // Integral -= Ex->GetBinContent(xbin+1,Y1+1,Z1+1)*Unit0;

                            if(abs(Ex->GetBinContent(xbin+1,Y0+1,Z0+1)) > 1E35 ){
                                cout<<"x: "<<xbin<<"; y: "<< Y0<<"; z: "<<Z0<<endl;
                                Unknown_jump = true;
                                // goto nextEntry;
                            }

                            if(abs(Ex->GetBinContent(xbin+1,Y1+1,Z1+1)) > 1E35){
                                cout<<"x: "<<xbin<<"; y: "<< Y1<<"; z: "<<Z1<<endl;
                                Unknown_jump = true;
                                // goto nextEntry;
                            }

                            if(abs((Ex->GetBinContent(xbin+1,Y0+1,Z0+1)-E0)*Unit0) > 0.5 ){
                                cout<<"xx big gap, x: "<<xbin<<"; y: "<< Y0<<"; z: "<<Z0<<"; gap: "<< (Ex->GetBinContent(xbin+1,Y0+1,Z0+1)-E0)*Unit0<<endl;
                            }

                            if(abs((Ex->GetBinContent(xbin+1,Y1+1,Z1+1)-E0)*Unit0) > 0.5){
                                cout<<"xx big gap, x: "<<xbin<<"; y: "<< Y1<<"; z: "<<Z1<<"; gap: "<< (Ex->GetBinContent(xbin+1,Y1+1,Z1+1)-E0)*Unit0<<endl;
                            }



                        }
                        if(Unknown_jump){continue;}

                        skipX:

                        if(Yend<Ymiddle){Y0 = Yend; Y1 = Ymiddle; X0 = Xend; X1 = Xmiddle; Z0 = Zend; Z1 = Zmiddle; }
                        if(Yend>Ymiddle){Y0 = Ymiddle; Y1 = Yend; X0 = Xmiddle; X1 = Xend; Z0 = Zmiddle; Z1 = Zend; }
                        if(Yend == Ymiddle){goto skipY;}

                        for (int ybin = Y0; ybin < Y1 + 1; ybin++){

                            if(ybin == Y0 || ybin == Y1){
                                Unit1 = DetectorReso[1]*0.5;
                            }
                            else{
                                Unit1 = DetectorReso[1];
                            }

                            Integral += Ey->GetBinContent(X1+1,ybin+1,Z0+1)*Unit1;
                            Integral -= Ey->GetBinContent(X0+1,ybin+1,Z1+1)*Unit1;


                            if(abs(Ey->GetBinContent(X1+1,ybin+1,Z0+1)) > 1E35  ){
                                cout<<"x: "<<X1<<"; y: "<< ybin<<"; z: "<<Z0<<endl;
                                Unknown_jump = true;
                                // goto nextEntry;
                            }

                            if(abs(Ey->GetBinContent(X0+1,ybin+1,Z1+1)) > 1E35){
                                cout<<"x: "<<X0<<"; y: "<< ybin<<"; z: "<<Z1<<endl;
                                Unknown_jump = true;
                                // goto nextEntry;
                            }

                            if(abs(Ey->GetBinContent(X1+1,ybin+1,Z0+1)*Unit1) > 0.5  ){
                                cout<<"yy big gap, x: "<<X1<<"; y: "<< ybin<<"; z: "<<Z0<<endl;
                            }

                            if(abs(Ey->GetBinContent(X0+1,ybin+1,Z1+1)*Unit1) > 0.5){
                                cout<<"yy big gap, x: "<<X0<<"; y: "<< ybin<<"; z: "<<Z1<<endl;
                            }


                        }

                        if(Unknown_jump){continue;}

                        skipY:

                        if(Zend<Zmiddle){Z0 = Yend; Z1 = Ymiddle; X0 = Xend; X1 = Xmiddle; Y0 = Yend; Y1 = Ymiddle; }
                        if(Zend>Zmiddle){Z0 = Ymiddle; Z1 = Yend; X0 = Xmiddle; X1 = Xend; Y0 = Ymiddle; Y1 = Yend; }
                        if(Zend == Zmiddle){goto skipZ;}

                        for (int zbin = Z0; zbin < Z1 + 1; zbin++){

                            if(zbin == Z0 || zbin == Z1){
                                Unit2 = DetectorReso[2]*0.5;
                            }
                            else{
                                Unit2 = DetectorReso[2];
                            }

                            Integral += Ez->GetBinContent(X1+1,Y1+1,zbin+1)*Unit2;
                            Integral -= Ez->GetBinContent(X0+1,Y0+1,zbin+1)*Unit2;

                            if(abs(Ez->GetBinContent(X1+1,Y1+1,zbin+1)) > 1E35){
                                cout<<"x: "<<X1<<"; y: "<< Y1<<"; z: "<<zbin<<endl;
                                Unknown_jump = true;
                            }

                            if(abs(Ez->GetBinContent(X0+1,Y0+1,zbin+1)) > 1E35){
                                cout<<"x: "<<X0<<"; y: "<< Y0<<"; z: "<<zbin<<endl;
                                Unknown_jump = true;
                            }

                            if(abs(Ez->GetBinContent(X1+1,Y1+1,zbin+1)*Unit2) > 0.5){
                                cout<<"zz big gap, x: "<<X1<<"; y: "<< Y1<<"; z: "<<zbin<<endl;
                            }

                            if(abs(Ez->GetBinContent(X0+1,Y0+1,zbin+1)*Unit2)> 0.5){
                                cout<<"zz big gap, x: "<<X0<<"; y: "<< Y0<<"; z: "<<zbin<<endl;
                            }

                        }

                    }
                    else
                    {
                        //what to do?????
                    }

                }
                if(abs(Efield[zbin + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][2]) > 0.5*float_max){
                    if(abs(Efield[(zbin+1) + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][2]) < 0.5*float_max||
                       abs(Efield[(zbin-1) + ybin * Resolution[2] + xbin * Resolution[2] * Resolution[1]][2]) < 0.5*float_max)
                    {
                        //hope
                        //maybe one unknown variable or 2.. safe to take two path
                    }
                    else
                    {
                        //what to do?????
                    }
                }
                else{
                    // Loop over all coordinates dx,dy,dz
                    for (unsigned coord = 0; coord < 3; coord++) {
                        // Fill interpolated grid points into histograms. bin=0 is underflow, bin = nbin+1 is overflow
                        Emap[coord].SetBinContent(xbin + 1, ybin + 1, zbin + 1, Efield[zbin + ybin * Resolution[2] +
                                                                                       xbin * Resolution[2] *
                                                                                       Resolution[1]][coord]);
                    } // end coordinate loop
                }
//                std::cout<<"xbin: "<<xbin<<"; ybin: "<<ybin<<"; zbin: "<<zbin<<"---Ex: "<<Efield[zbin+ybin*Resolution[2]+xbin*Resolution[2]*Resolution[1]][0]<<"; Ey: "<<Efield[zbin+ybin*Resolution[2]+xbin*Resolution[2]*Resolution[1]][1]<<"; Ez: "<< Efield[zbin+ybin*Resolution[2]+xbin*Resolution[2]*Resolution[1]][2]<<std::endl;
            } // end zbin loop
        } // end ybin loop
    } // end zbin loop

}

*/


