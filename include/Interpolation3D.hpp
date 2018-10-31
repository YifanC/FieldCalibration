// #define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <tuple>
#include <array>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <string>
#include <iomanip>
#include <chrono>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include "ThreeVector.hpp"
#include "LaserTrack.hpp"
#include "Matrix3x3.hpp"

#ifndef INTERPOLATION3D_H
#define INTERPOLATION3D_H

typedef CGAL::Exact_predicates_inexact_constructions_kernel InterpKernel;
typedef CGAL::Triangulation_3<InterpKernel> Triangulation;
typedef CGAL::Triangulation_vertex_base_with_info_3<std::pair<unsigned long, unsigned long>, InterpKernel> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
// Use the Fast_location tag. Default or Compact_location works too.
// typedef CGAL::Delaunay_triangulation_3<InterpKernel, Tds, CGAL::Fast_location > Delaunay;
// Use the Lock for parallel computing
typedef CGAL::Delaunay_triangulation_3<InterpKernel, Tds> Delaunay;
typedef Delaunay::Point Point;

std::vector<std::pair<unsigned, float>> GetClosestTracksInfo(std::vector<LaserTrack> &, const unsigned);

std::vector<std::pair<unsigned int, unsigned int>> GetClosestLaserSample(std::vector<LaserTrack> &, const unsigned);

std::array<float, 2> AnglesFromPoynting(ThreeVector<float> &);

ThreeVector<float> PoyntingFromAngles(const std::array<float, 2> &);

bool PairSortFunction(std::pair<unsigned int, float>, std::pair<unsigned int, float>);

Delaunay TrackMesher(const std::vector<LaserTrack> &);

Point VectorToPoint(ThreeVector<float> &);

ThreeVector<float> PointToVector(Point &);

//ThreeVector<float>
//InterpolateCGAL(const std::vector<LaserTrack> &LaserTrackSet, const std::vector<LaserTrack> &LaserMeshSet,
//                const Delaunay &Mesh, ThreeVector<float> Location);
//
//std::vector<ThreeVector<float>>
//InterpolateMap(const std::vector<LaserTrack> &LaserTrackSet, const std::vector<LaserTrack> &LaserMeshSet,
//               const Delaunay &Mesh, const TPCVolumeHandler &TPC);

ThreeVector<float>
InterpolateCGAL(const std::vector<LaserTrack> &LaserTrackSet, const Delaunay &Mesh, ThreeVector<float> Location);

std::vector<ThreeVector<float>>
InterpolateMap(const std::vector<LaserTrack> &LaserTrackSet, const Delaunay &Mesh, const TPCVolumeHandler &TPC);

void InterpolateTrack(LaserTrack &, const std::vector<LaserTrack> &, const Delaunay &);

/////////////////E map session
typedef CGAL::Triangulation_vertex_base_with_info_3<int, InterpKernel> xVb;
typedef CGAL::Triangulation_data_structure_3<xVb> xTds;
typedef CGAL::Delaunay_triangulation_3<InterpKernel, xTds> xDelaunay;
typedef xDelaunay::Point xPoint;

xDelaunay Mesher(std::vector<ThreeVector<float>> &Position, TPCVolumeHandler &TPC);

xPoint xVectorToPoint(ThreeVector<float> &);

std::pair<ThreeVector<float>, ThreeVector<float>>
EInterpolateCGAL(std::vector<ThreeVector<float>> &vn, std::vector<ThreeVector<float>> &En, std::vector<ThreeVector<float>> &Position, const xDelaunay &Mesh,
                 ThreeVector<float> Location, const TPCVolumeHandler &TPC);

float EcompInterpolateCGAL(std::vector<float> &En, std::vector<ThreeVector<float>> &Position,
                           const xDelaunay &Mesh, ThreeVector<float> Location, const TPCVolumeHandler &TPC);

std::pair<std::vector<ThreeVector<float>>, std::vector<ThreeVector<float>>>
EInterpolateMap(std::vector<ThreeVector<float>> &vn, std::vector<ThreeVector<float>> &En, std::vector<ThreeVector<float>> &Position, const xDelaunay &Mesh,
                const TPCVolumeHandler &TPC, ThreeVector<unsigned long> EReso);

std::vector<ThreeVector<float>>
EcompInterpolateMap(std::vector<float> &Ex, std::vector<ThreeVector<float>> &PositionX, const xDelaunay &MeshX,
                    std::vector<float> &Ey, std::vector<ThreeVector<float>> &PositionY, const xDelaunay &MeshY,
                    std::vector<float> &Ez, std::vector<ThreeVector<float>> &PositionZ, const xDelaunay &MeshZ,
                    const TPCVolumeHandler &TPC, ThreeVector<unsigned long> EReso);




/////////////////E map session

#endif