#include "../include/LaserTrack.hpp"

// Default constructor
LaserTrack::LaserTrack() {

}

// Constructor used for reconstructed laser data
LaserTrack::LaserTrack(const TVector3 &InEntryPoint, const TVector3 &InExitPoint,
                       const std::vector<TVector3> &RecoTrack) {
    // Resever enogh space for the LaserReco vector. Can be important if you go to the Memory limit
    LaserReco.reserve(RecoTrack.size());

    // Store entry and exit point
    EntryPoint = ThreeVector<float>(InEntryPoint);
    ExitPoint = ThreeVector<float>(InExitPoint);

    // Loop over all RecoTrack points
    for (const auto &TrackPoint : RecoTrack) {
        // Fill TrackPoints as vectors into the LaserReco vector
        LaserReco.push_back(ThreeVector<float>(TrackPoint));
    }

    // Initialize LaserDisplacement vector
    LaserDisplacement.resize(LaserReco.size());
}

// Constructor used for reconstructed laser data
LaserTrack::LaserTrack(const ThreeVector<float> &InEntryPoint, const ThreeVector<float> &InExitPoint,
                       const std::vector<ThreeVector<float>> &RecoTrack) {
    // Resever enogh space for the LaserReco vector. Can be important if you go to the Memory limit
    LaserReco.reserve(RecoTrack.size());

    // Store entry and exit point
    EntryPoint = InEntryPoint;
    ExitPoint = InExitPoint;

    // Loop over all RecoTrack points
    for (const auto &TrackPoint : RecoTrack) {
        // Fill TrackPoints as vectors into the LaserReco vector
        LaserReco.push_back(TrackPoint);
    }

    // Initialize LaserDisplacement vector
    LaserDisplacement.resize(LaserReco.size());
}

LaserTrack::LaserTrack(std::array<float, 2> &Angles, ThreeVector<float> &Position, const TPCVolumeHandler &TPCVolume)
        : TrackAngles(Angles), LaserPosition(Position) {
    FindBoundaries(TPCVolume);
}

////////////////////////////////////
//TODO:Check if it's right
// This part is to disguise the anode zero distortion as tracks, so they can be easily merged after iteration calculation
LaserTrack::LaserTrack(const std::vector<ThreeVector<float>> &RecoTrack, const std::vector<ThreeVector<float>> &LaserDisp) {

    if(RecoTrack.size()!=LaserDisp.size()){
        std::cerr << "ERROR: Please provide RecoTrack and LaserDisp in same size !" << std::endl;
    }

    // Resever enogh space for the LaserReco vector. Can be important if you go to the Memory limit
    LaserReco.reserve(RecoTrack.size());

    // Loop over all RecoTrack points
    for (const auto &TrackPoint : RecoTrack) {
        // Fill TrackPoints as vectors into the LaserReco vector
        LaserReco.push_back(TrackPoint);
    }

    // Initialize LaserDisplacement vector
    LaserDisplacement.resize(LaserReco.size());

    // Loop over all RecoTrack points
    for (const auto &DispVec : LaserDisp) {
        // Fill TrackPoints as vectors into the LaserReco vector
        LaserDisplacement.push_back(DispVec);
    }
}


//LaserTrack::LaserTrack(std::array<float, 2> &Angles, ThreeVector<float> &Position, const TPCVolumeHandler &TPCVolume)
//        : TrackAngles(Angles), LaserPosition(Position) {
//    FindBoundaries(TPCVolume);
//}

// This Constructor initializes the laser track vectors and fills them according to the angles and segments
LaserTrack::LaserTrack(const unsigned int NSegments, std::array<float, 2> &Angles, ThreeVector<float> &Position,
                       const TPCVolumeHandler &TPCVolume) : NumberOfTracksegments(NSegments), TrackAngles(Angles),
                                                            LaserPosition(Position) {
//   LaserTrue.resize(NumberOfTracksegments+1);
    LaserReco.resize(NumberOfTracksegments + 1);
    LaserDisplacement.resize(NumberOfTracksegments + 1);

    FindBoundaries(TPCVolume);
    FillTrack();
}

ThreeVector<float> LaserTrack::GetPoyntingVector() {
    // Create unit vector of Poynting vector from spherical coordinates
    ThreeVector<float> vec_res;
    vec_res[0] = (float) (std::cos(TrackAngles[0]) * std::sin(TrackAngles[1])); // cos(theta)*cos(phi)
    vec_res[1] = (float) std::sin(TrackAngles[0]); // sin(theta)
    vec_res[2] = (float) (std::cos(TrackAngles[0]) * std::cos(TrackAngles[1])); // cos(theta)*cos(phi)

    return vec_res;
}

void LaserTrack::FillTrack() {
    ThreeVector<float> Poynting = GetPoyntingVector();

    LaserReco.front() = EntryPoint;
    LaserReco.back() = ExitPoint;

    // Parameter length of a single track segment
    float SegmentParameter = ThreeVector<float>::VectorNorm(ExitPoint - EntryPoint) / (float) NumberOfTracksegments;

    // Fill segments
    for (unsigned segment = 1; segment < NumberOfTracksegments; segment++) {
        LaserReco[segment] = EntryPoint + (SegmentParameter * segment) * Poynting;
    }
}

void LaserTrack::FindBoundaries(const TPCVolumeHandler &TPCVolume) {
    ThreeVector<float> Poynting = GetPoyntingVector();

    // This is the absolute error allowed for the boundary conditions of the laser tracks
    float epsilon_abs = 1e-6;

//   TPCVolumeHandler TPCVolume(DetectorSize, DetectorOffset);
    auto SurfaceNormVec = TPCVolume.GetNormalVectors();
    auto NormVecOffset = TPCVolume.GetNormVectorOffset();

    std::vector<float> BoundaryParameter;
    std::vector<ThreeVector<float>> BoundaryPoint;

    for (unsigned NSurface = 0; NSurface < SurfaceNormVec.size(); NSurface++) {
        float IntersecParameter;

        // Calculate the intersection point (UnitVector is added and suptracted to prevent rounding errors)
        ThreeVector<float> IntersecPoint = LaserPosition + Poynting * IntersecParameter;
        bool IntersectionFlag = true;

        for (unsigned coord = 0; coord < IntersecPoint.size(); coord++) {
            // Check if Intersection Point is within the surface area with an error of epsilon_abs
            if (IntersecPoint[coord] >
                TPCVolume.GetDetectorSize()[coord] + TPCVolume.GetDetectorOffset()[coord] + epsilon_abs ||
                IntersecPoint[coord] < TPCVolume.GetDetectorOffset()[coord] - epsilon_abs || IntersecParameter < 0.0) {
                IntersectionFlag = false;
                break;
            }
        }
        if (IntersectionFlag) {
            BoundaryParameter.push_back(IntersecParameter);
            BoundaryPoint.push_back(IntersecPoint);
        }
    }

    // If there are no hits on the TPC the laser track will be initialized empty with an error
    if (!BoundaryParameter.size()) {
//     LaserTrue.clear();
//     LaserReco.clear();
        std::cout << "ERROR: Laser beam with start position = (" << LaserPosition[0] << ", " << LaserPosition[1] << ", "
                  << LaserPosition[2] << ") "
                  << "and track angles (theta, phi) = (" << TrackAngles[0] / M_PI * 180.0 << ", "
                  << TrackAngles[1] / M_PI * 180.0 << ") "
                  << "is missing the chamber!" << std::endl;
        EntryPoint = {0.0, 0.0, 0.0};
        ExitPoint = {0.0, 0.0, 0.0};
        return;
    } else if (BoundaryParameter.size() <
               2) // If the laser beam is produced in the volume add the laser position and the parameter as boundary conditions
    {
        BoundaryParameter.push_back(0.0);
        BoundaryPoint.push_back(LaserPosition);
        std::cout << "Warning: Laser beam with start position = (" << LaserPosition[0] << ", " << LaserPosition[1]
                  << ", " << LaserPosition[2] << ") "
                  << "and track angles (theta, phi) = (" << TrackAngles[0] / M_PI * 180.0 << ", "
                  << TrackAngles[1] / M_PI * 180.0 << ") "
                  << "is placed in the chamber!" << std::endl;
    }

    // Check order of boundary points and swap if they are in wrong order
    if (BoundaryParameter[0] > BoundaryParameter[1]) {
        std::swap(BoundaryPoint.front(), BoundaryPoint.back());
        std::swap(BoundaryParameter.front(), BoundaryParameter.back());
    }

    // Fill the entry and exit point
    EntryPoint = BoundaryPoint.front();
    ExitPoint = BoundaryPoint.back();
}

// This is the displacement algorithm using the track derivative method
void LaserTrack::DerivativeDisplAlgo() {
    // Initialize scaling parameter for point on true laser beam
    float TrueTrackPara;

    // Vector which points from entry to exit point
    ThreeVector<float> TrueTrack = ExitPoint - EntryPoint;

    // Only Displacement calculation for intermediate points if there are 3 or more track points
    if (LaserReco.size() > 2) {
        // Calculate displacement of all intermediate track samples with the derivative method
        for (unsigned sample_no = 1; sample_no < LaserReco.size() - 1; sample_no++) {
            // Calculate intersection parameter on the true track (Parameter*TrueTrackVector = IntersectionPoint)
            TrueTrackPara = ThreeVector<float>::DotProduct(LaserReco[sample_no - 1] - LaserReco[sample_no + 1],
                                                           LaserReco[sample_no] - EntryPoint)
                            / ThreeVector<float>::DotProduct(LaserReco[sample_no - 1] - LaserReco[sample_no + 1],
                                                             TrueTrack);

            // Calculate vector between IntersectionPoint and LaserRecoPoint, which is the displacement vector
            LaserDisplacement[sample_no] = LaserReco[sample_no] - (EntryPoint + TrueTrackPara * TrueTrack);
        }
    } // end if

    // If there are two and more points, shift entry and exit point to true position 
    if (LaserReco.size() > 1) {
        // Displacement for first and last track samples
        LaserDisplacement.front() = LaserReco.front() - EntryPoint;
        LaserDisplacement.back() = LaserReco.back() - ExitPoint;
    }
}

// This algorithm calculates the displacement by projecting the reco samples to the closest point on the true track
//void LaserTrack::ClosestPointDisplAlgo(bool CorrMapFlag = true, int Nstep = 1) {
void LaserTrack::ClosestPointDisplAlgo(int Nstep = 1) {
    // Initialize scaling parameter for point on true laser beam
    float TrueTrackPara;

    // Vector which points from entry to exit point
    ThreeVector<float> TrueTrack = ExitPoint - EntryPoint;

    for (unsigned sample_no = 0; sample_no < LaserReco.size(); sample_no++) {
        // Calculate parameter corresponding to the closest distance between true track and reco track sample
        TrueTrackPara = ThreeVector<float>::DotProduct(LaserReco[sample_no] - EntryPoint, TrueTrack) /
                        ThreeVector<float>::DotProduct(TrueTrack, TrueTrack);

        // Calculate displacement
//        if (CorrMapFlag) {
            LaserDisplacement[sample_no] = (EntryPoint + TrueTrackPara * TrueTrack) - LaserReco[sample_no];
            LaserDisplacement[sample_no] /= Nstep;
//        } else {
//            LaserDisplacement[sample_no] = LaserReco[sample_no] - (EntryPoint + TrueTrackPara * TrueTrack);
//        }
    }
}

// Warning: There's some problem with this algorithm...
// This algorithm uses first the closest point displacement algorithm. Then it shifts the reco points linearly along the true track to equalize the reco and truth track length 
void LaserTrack::LinearStretchDisplAlgo(bool CorrMapFlag = true) {
    // Execute closest point algorithm
    ClosestPointDisplAlgo();

    // Vector which points from entry to exit point
    ThreeVector<float> TrueTrack = ExitPoint - EntryPoint;

    // Stretch vector for every sample
    std::vector<ThreeVector<float>> Stretch;
    Stretch.reserve(GetNumberOfSamples());

    // Initialize vector with corrected points
    std::vector<ThreeVector<float>> Corrected = LaserReco;

    // Loop over all samples of corrected
    for (unsigned long sample = 0; sample < Corrected.size(); sample++) {
        Corrected[sample] -= LaserDisplacement[sample];
    }

    // Calculate EntryPoint and ExitPoint shift
    float EntryPointShift = ThreeVector<float>::DotProduct(TrueTrack.GetUnitVector(), EntryPoint - Corrected.front());
    float ExitPointShift = ThreeVector<float>::DotProduct(TrueTrack.GetUnitVector(), ExitPoint - Corrected.back());

    // Loop over all track samples (needs to be seperate, in order to not interfere with the stretch result)
    for (unsigned long sample = 0; sample < Corrected.size(); sample++) {
        // Calculate the stretch for every sample point
        Stretch.push_back((ExitPointShift - EntryPointShift) /
                          ThreeVector<float>::VectorNorm(Corrected.back() - Corrected.front()) *
                          (Corrected.front() - Corrected[sample])
                          - EntryPointShift * TrueTrack.GetUnitVector());
    }

    // Loop over all track samples
    for (unsigned long sample = 0; sample < GetNumberOfSamples(); sample++) {
        // Add stretch correction (subtract stretch displacement) to existing displacement
        LaserDisplacement[sample] -= Stretch[sample];
    }

}


void LaserTrack::DistortTrack(std::string MapFileName, const TPCVolumeHandler &TPCVolume) {
    float epsilon_abs = 1e-6;

    TFile *FieldFile = new TFile(MapFileName.c_str());
    if (FieldFile->IsZombie()) {
        std::cout << "ERROR: " << MapFileName << " can not be opened!" << std::endl;
        return;
    }

    std::vector<TH3F *> DistortMap;
    DistortMap.push_back((TH3F *) FieldFile->Get("Distortion_Field_X"));
    DistortMap.push_back((TH3F *) FieldFile->Get("Distortion_Field_Y"));
    DistortMap.push_back((TH3F *) FieldFile->Get("Distortion_Field_Z"));

    float InterpolDistortion = 0.0;

    // Add field distortion to every track segment point
    for (unsigned segment = 0; segment < NumberOfTracksegments + 1; segment++) {
        for (unsigned coord = 0; coord < LaserReco[segment].size(); coord++) {
            if (InterpolDistortion == DistortMap[coord]->Interpolate(LaserReco[segment][0], LaserReco[segment][1] +
                                                                                            TPCVolume.GetDetectorSize()[1] /
                                                                                            (float) 2.0,
                                                                     LaserReco[segment][2])) {
                LaserReco[segment][coord] += InterpolDistortion;
            }
        }
    }
    for (int coord = 0; coord < LaserReco[0].size(); coord++) delete DistortMap[coord];
    FieldFile->Close();
    delete FieldFile;
    gDirectory->GetList()->Delete();

//   for(unsigned segment = 0; segment < NumberOfTracksegments+1; segment++)
//   {
//     std::cout << LaserTrue[segment][0] << " " << LaserTrue[segment][1] << " " << LaserTrue[segment][2] << std::endl;
//     std::cout << LaserReco[segment][0] << " " << LaserReco[segment][1] << " " << LaserReco[segment][2] << std::endl;
//   }
}

// This function calculates the track displacement
void LaserTrack::CalcDisplacement(const DisplacementAlgo &Algo, int Nstep) {
    // Switch to chosen displacement algorithm.
    // In case you add algorithms introduce more cases
    switch (Algo) {
        case TrackDerivative :
            DerivativeDisplAlgo();
            break;
//        case ClosestPointCorr :
//            ClosestPointDisplAlgo(true, Nstep);
//            break;
//        case ClosestPointDist :
//            ClosestPointDisplAlgo(false);
//            break;
        case ClosestPoint :
            ClosestPointDisplAlgo(Nstep);
            break;
        case LinearStretch :
            LinearStretchDisplAlgo();
            break;
        default :
            ClosestPointDisplAlgo(Nstep);
            break;
    }
}

//// Add correction to the reco position (correction is the negative displacement).
//// This is important for generating a displacement map in non-distorted detector coordinates.
//void LaserTrack::AddCorrectionToRecoM() {
//    // Loop over track points
//    for (unsigned long sample = 0; sample < LaserReco.size(); sample++) {
//        // "Add" DISTORTION to reconstructed track position
//        // The sign should be minus here, since AddCorrectionToReco is only used for meshing on true coord
//        // In that case, we calculate the distortion vectors
//        // To drag the reco track back on true, we use minus sign for the direction of the vectors
//        LaserReco[sample] -= LaserDisplacement[sample];
//    }
//}
//
//// Add correction to the reco position (correction is the negative displacement).
//// This is important for generating a displacement map in non-distorted detector coordinates.
//void LaserTrack::AddCorrectionToRecoP() {
//    // Loop over track points
//    for (unsigned long sample = 0; sample < LaserReco.size(); sample++) {
//        // Add CORRECTION to reconstructed track position
//        // drag to the reco to the true postion
//        LaserReco[sample] += LaserDisplacement[sample];
//    }
//}

// Add correction to the reco position (correction is the negative displacement).
// This is important for generating a displacement map in non-distorted detector coordinates.
void LaserTrack::AddCorrectionToReco() {
    // Loop over track points
    for (unsigned long sample = 0; sample < LaserReco.size(); sample++) {
        // Add CORRECTION to reconstructed track position
        // drag to the reco to the true postion
        LaserReco[sample] += LaserDisplacement[sample];
    }
}

// Add correction to the reco position with (e.g 1/10 of) the correction vector calculated by the other side laser displacement
void LaserTrack::AddCorrectionToRecoPart(std::vector<ThreeVector<float>> CorrPart) {
    if (LaserReco.size() == CorrPart.size()) {
        for (unsigned long sample = 0; sample < LaserReco.size(); sample++) {
            // Add part of the correction to reconstructed track position
            LaserReco[sample] += CorrPart[sample];
        }
    }
    else {
        std::cerr << "This laser track does not have the same size as the corresponding correction vector. Please check." << std::endl;
    }
}

// Calculating the displacement from the true track to the reco track
// Since at this stage, the laser tracks are dragged to the true track lines already, the LaserReco is on true track lines
// Run it with the lasersets and supply the origin reco tracks
void LaserTrack::Displacement(LaserTrack LaserTrackReco, bool Corr = true) {
    // Check if the tracks of the two LaserSets have the same size
    if (LaserReco.size() == LaserTrackReco.LaserReco.size()) {
        // LaserReco is the true laser tracks
        if (Corr) {
            // Loop over track points
            for (unsigned long sample = 0; sample < LaserReco.size(); sample++) {
                // Calculate the correction vector
                LaserDisplacement[sample] = LaserReco[sample] - LaserTrackReco.LaserReco[sample];
            }
        } else {
            for (unsigned long sample = 0; sample < LaserReco.size(); sample++) {
                // Calculate the distortion vector
                LaserDisplacement[sample] = LaserTrackReco.LaserReco[sample] - LaserReco[sample];
            }
        }
    } else {
        std::cerr << "ERROR: Please supply the same size of laser sets to calculate the displacement" << std::endl;
    }

}


void LaserTrack::AddToDisplacement(ThreeVector<float> &AdditionVector, unsigned long sample) {
    LaserDisplacement[sample] += AdditionVector;
}

ThreeVector<float> LaserTrack::GetLaserPosition() {
    return LaserPosition;
}

std::array<float, 2> LaserTrack::GetAngles() {
    return TrackAngles;
}

unsigned long LaserTrack::GetNumberOfSamples() const {
    return LaserReco.size();
}

ThreeVector<float> LaserTrack::GetSamplePosition(const unsigned int &SampleNumber) const {
    return LaserReco[SampleNumber];
}

ThreeVector<float> LaserTrack::GetDisplacement(const unsigned int &SampleNumber) const {
    return LaserDisplacement[SampleNumber];
}

ThreeVector<float> LaserTrack::GetEntryPoint() {
    return EntryPoint;
}

ThreeVector<float> LaserTrack::GetExitPoint() {
    return ExitPoint;
}

void LaserTrack::AppendSample(ThreeVector<float> &SamplePosition, ThreeVector<float> &SampleDisplacement) {
    LaserReco.push_back(SamplePosition);
    LaserDisplacement.push_back(SampleDisplacement);
}

void LaserTrack::AppendSample(ThreeVector<float> &SamplePosition) {
    LaserReco.push_back(SamplePosition);
    LaserDisplacement.push_back(ThreeVector<float>(0.0, 0.0, 0.0));
}

void LaserTrack::AppendSample(float SamplePos_x, float SamplePos_y, float SamplePos_z) {
    LaserReco.push_back(ThreeVector<float>(SamplePos_x, SamplePos_y, SamplePos_z));
    LaserDisplacement.push_back(ThreeVector<float>(0.0, 0.0, 0.0));
}

void LaserTrack::AppendSample(float SamplePos_x, float SamplePos_y, float SamplePos_z, float SampleDispl_x,
                              float SampleDispl_y, float SampleDispl_z) {
    LaserReco.push_back(ThreeVector<float>(SamplePos_x, SamplePos_y, SamplePos_z));
    LaserDisplacement.push_back(ThreeVector<float>(SampleDispl_x, SampleDispl_y, SampleDispl_z));
}

void LaserTrack::DistortTracks(std::vector<LaserTrack> &LaserTracks, const std::string &MapFileName,
                               const TPCVolumeHandler &TPCVolume) {
    float epsilon_abs = 1e-6;

    TFile *FieldFile = new TFile(MapFileName.c_str());
    if (FieldFile->IsZombie()) {
        std::cout << "ERROR: " << MapFileName << " can not be opened!" << std::endl;
        return;
    }

    std::vector<TH3F *> DistortMap;
    DistortMap.push_back((TH3F *) FieldFile->Get("Distortion_Field_X"));
    DistortMap.push_back((TH3F *) FieldFile->Get("Distortion_Field_Z"));
    DistortMap.push_back((TH3F *) FieldFile->Get("Distortion_Field_Y"));

    float InterpolDistortion = 0.0;

    for (unsigned track_no = 0; track_no < LaserTracks.size(); track_no++) {
//     // Add field distortion to every track segment point
        for (unsigned segment = 0; segment < LaserTracks[track_no].NumberOfTracksegments + 1; segment++) {
            for (unsigned coord = 0; coord < LaserTracks[track_no].LaserReco[segment].size(); coord++) {
                LaserTracks[track_no].LaserReco[segment][coord] += DistortMap[coord]->Interpolate(
                        LaserTracks[track_no].LaserReco[segment][0] - TPCVolume.GetDetectorOffset().at(0),
                        LaserTracks[track_no].LaserReco[segment][1] - TPCVolume.GetDetectorOffset().at(1),
                        LaserTracks[track_no].LaserReco[segment][2] - TPCVolume.GetDetectorOffset().at(2)
                );
// 	  std::cout << LaserTracks[track_no].LaserTrue[segment][coord] << " " << LaserTracks[track_no].LaserReco[segment][coord] << " ";
// 	  std::cout << DistortMap[coord]->Interpolate(LaserTracks[track_no].LaserTrue[segment][0],LaserTracks[track_no].LaserTrue[segment][1]+DetectorSize[1]/(float)2.0,LaserTracks[track_no].LaserTrue[segment][2]) << " ";
            }
        }
//       std::cout << std::endl;
    }
//   
    for (unsigned coord = 0; coord < DistortMap.size(); coord++) {
        delete DistortMap[coord];
    }
    FieldFile->Close();
    delete FieldFile;
    gDirectory->GetList()->Delete();
}

std::vector<ThreeVector<float>> LaserTrack::GetReco() {
    return LaserReco;
}

ThreeVector<float> LaserTrack::GetFront() const {
    return LaserReco.front();
}

ThreeVector<float> LaserTrack::GetBack() const {
    return LaserReco.back();
}
