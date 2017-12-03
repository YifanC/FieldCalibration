#include "../include/Laser.hpp"

Laser::Laser() {}

Laser::Laser(std::vector<LaserTrack> InputSet) {
    LaserTrackSet = InputSet;
}

std::vector<LaserTrack>::iterator Laser::begin() {
    return LaserTrackSet.begin();
}

std::vector<LaserTrack>::iterator Laser::end() {
    return LaserTrackSet.end();
}

void Laser::AppendTrack(const LaserTrack &InputTrack) {
    LaserTrackSet.push_back(InputTrack);
}

// Merges a vector of Lasers into one Laser object
Laser Laser::Merge(std::vector<Laser> &LaserVec) {
    // Initialize output of type Laser
    Laser MergedLaserSets;

    // Loop over all input vector entries
    for (unsigned long index = 0; index < LaserVec.size(); index++) {
        MergedLaserSets.LaserTrackSet.insert(MergedLaserSets.LaserTrackSet.end(),
                                             LaserVec.at(index).LaserTrackSet.begin(),
                                             LaserVec.at(index).LaserTrackSet.end());
    }

    LaserVec.clear();
}

// Merge LaserB into current Laser and rewrite the current Laser
void Laser::Mergewith(const Laser &LaserB) {
    std::vector<LaserTrack> LaserVec = LaserB.GetTrackSet();
    for (unsigned long n = 0; n < LaserVec.size(); n++){
        LaserTrackSet.push_back(LaserVec[n]);
    }
}


LaserTrack Laser::GetTrack(const long unsigned int &TrackNumber) {
    LaserTrackSet.at(TrackNumber);
}

LaserTrack Laser::GetLastTrack() {
    return LaserTrackSet.back();
}

LaserTrack Laser::GetFirstTrack() {
    return LaserTrackSet.front();
}

long unsigned int Laser::GetNumberOfTracks() const {
    return LaserTrackSet.size();
}

std::vector<LaserTrack> Laser::GetTrackSet() const {
    return LaserTrackSet;
}

void Laser::DistortTrackSet(std::string MapFileName, TPCVolumeHandler &TPCVolume) {
    LaserTrack::DistortTracks(LaserTrackSet, MapFileName, TPCVolume);
}

// Applies correction algorithm on all tracks of a laser track set
void Laser::CalcDisplacement(const LaserTrack::DisplacementAlgo &Algo, int Nstep) {
    // Loop over all tracks, and calculate displacement
    for (auto &Track : LaserTrackSet) {
        Track.CalcDisplacement(Algo, Nstep);
    }
}

//// Add displacement to the reco position. This is important for generating a displacement map in non-distorted detector coordinates
//void Laser::AddCorrectionToReco(bool plus = true) {
//    // Loop over all tracks, and calculate displacement
//    for (auto &Track : LaserTrackSet) {
//        if (plus) { Track.AddCorrectionToRecoP(); }
//        if (!plus) { Track.AddCorrectionToRecoM(); }
//    }
//}

// Add displacement to the reco position. This is important for generating a displacement map in non-distorted detector coordinates
void Laser::AddCorrectionToReco() {
    // Loop over all tracks, and calculate displacement
    for (auto &Track : LaserTrackSet) {
        Track.AddCorrectionToReco();
    }
}

// Add displacement to the reco position. This is important for generating a displacement map in non-distorted detector coordinates
void Laser::SetDisplacement(Laser LaserReco, bool Corr = true) {
    // Loop over all tracks, and calculate displacement
    if (LaserTrackSet.size() == LaserReco.GetTrackSet().size()) {
        for (int i = 0; i < LaserTrackSet.size(); i++) {
            LaserTrackSet[i].Displacement(LaserReco.GetTrackSet()[i], Corr);
        }
    } else {
        std::cerr << "ERROR: Please supply the same size of laser sets to calculate the displacement" << std::endl;
    }

}

void Laser::InterpolateTrackSet(const std::vector<LaserTrack> &LaserTracks, const Delaunay &Mesh) {
    for (auto &Track : LaserTrackSet) {
        InterpolateTrack(Track, LaserTracks, Mesh);
    }
}

void Laser::InterpolateTrackSet(const Laser &LaserTracks, const Delaunay &Mesh) {
    InterpolateTrackSet(LaserTracks.GetTrackSet(), Mesh);
}

void Laser::DrawTrack(const long unsigned int &TrackNumber) {

    TApplication theApp("App", 0, 0);

//   std::cout << LaserTrackSet.at(2).GetNumberOfSamples() << std::endl;
    TGraph2D *TrueTrack = new TGraph2D(2);
    TPolyLine3D *DistortedTrack = new TPolyLine3D(LaserTrackSet.at(TrackNumber).GetNumberOfSamples());
    TPolyLine3D *CorrectedTrack = new TPolyLine3D(LaserTrackSet.at(TrackNumber).GetNumberOfSamples());


    for (unsigned sample_no = 0; sample_no < LaserTrackSet.at(TrackNumber).GetNumberOfSamples(); sample_no++) {
        DistortedTrack->SetPoint(sample_no,
                                 LaserTrackSet.at(TrackNumber).GetSamplePosition(sample_no).at(0),
                                 LaserTrackSet.at(TrackNumber).GetSamplePosition(sample_no).at(1),
                                 LaserTrackSet.at(TrackNumber).GetSamplePosition(sample_no).at(2));
        CorrectedTrack->SetPoint(sample_no,
                                 LaserTrackSet.at(TrackNumber).GetSamplePosition(sample_no).at(0) +
                                 LaserTrackSet.at(TrackNumber).GetDisplacement(sample_no).at(0),
                                 LaserTrackSet.at(TrackNumber).GetSamplePosition(sample_no).at(1) +
                                 LaserTrackSet.at(TrackNumber).GetDisplacement(sample_no).at(1),
                                 LaserTrackSet.at(TrackNumber).GetSamplePosition(sample_no).at(2) +
                                 LaserTrackSet.at(TrackNumber).GetDisplacement(sample_no).at(2));
    }

    TrueTrack->SetPoint(0, LaserTrackSet.at(TrackNumber).GetEntryPoint().at(0),
                        LaserTrackSet.at(TrackNumber).GetEntryPoint().at(1),
                        LaserTrackSet.at(TrackNumber).GetEntryPoint().at(2));
    TrueTrack->SetPoint(1, LaserTrackSet.at(TrackNumber).GetExitPoint().at(0),
                        LaserTrackSet.at(TrackNumber).GetExitPoint().at(1),
                        LaserTrackSet.at(TrackNumber).GetExitPoint().at(2));


    std::string PictureName = "Track_" + std::to_string(TrackNumber);

    TrueTrack->SetLineColor(kBlue);
    DistortedTrack->SetLineColor(kRed);
    CorrectedTrack->SetLineColor(kGreen);

    TCanvas *C0 = new TCanvas(PictureName.c_str(), PictureName.c_str(), 1000, 700);
    TrueTrack->Draw("LINE");
    DistortedTrack->Draw("same");
    CorrectedTrack->Draw("same");
    C0->Print((PictureName + ".png").c_str(), "png");

    theApp.Run();

    delete DistortedTrack;
    delete CorrectedTrack;
    delete TrueTrack;
    delete C0;
    gDirectory->GetList()->Delete();
}



