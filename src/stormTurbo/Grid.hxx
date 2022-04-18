/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person
/// obtaining a copy of this software and associated documentation
/// files (the "Software"), to deal in the Software without
/// restriction, including without limitation the rights  to use,
/// copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the
/// Software is furnished to do so, subject to the following
/// conditions:
///
/// The above copyright notice and this permission notice shall be
/// included in all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
/// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
/// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
/// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
/// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
/// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
/// OTHER DEALINGS IN THE SOFTWARE.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///

#pragma once

#include <iostream>
#include <functional>

#include <stormBase.hxx>
#include <stormSolvers/Mat.hxx>

#include <stormTurbo/Array.hxx>

namespace Storm::Turbo {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Structured 3D grid. 
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
class Grid3D final : public NonCopyable {
private:
  size_t StrokeWidth_;
  size_t NumCellsX_, NumCellsY_, NumCellsZ_;
  Array3D<Vec3D<real_t>> NodePos_;
  Array3D<real_t> CellVolume_;
  Array3D<Vec3D<real_t>> CellCenter_;
  Array3D<real_t> FaceAreaX_, FaceAreaY_, FaceAreaZ_;
  Array3D<Vec3D<real_t>> FaceNormalX_, FaceNormalY_, FaceNormalZ_;
  Array3D<Vec3D<real_t>> FaceCenterX_, FaceCenterY_, FaceCenterZ_;

public:

  /// @brief Initialize an empty 3D grid.
  Grid3D() : StrokeWidth_(0), NumCellsX_(0), NumCellsY_(0), NumCellsZ_(0) {}

  /// @brief Build a 3D grid with the specified size and coordinate transform.
  void Build(size_t numGhostCells,
             size_t numCellsX, size_t numCellsY, size_t numCellsZ,
             std::function<Vec3D<real_t>(Vec3D<size_t>)> const& nodeTransform = {});

  /// @brief Total number of cells in the X direction.
  size_t SizeX() const noexcept {
    return NumCellsX_ + 2*StrokeWidth_;
  }

  /// @brief Total number of cells in the Y direction.
  size_t SizeY() const noexcept {
    return NumCellsY_ + 2*StrokeWidth_;
  }
  
  /// @brief Total number of cells in the Z direction.
  size_t SizeZ() const noexcept {
    return NumCellsZ_ + 2*StrokeWidth_;
  }

  /// @brief Get the first X-index of the interior cells. 
  size_t BeginCellX() const noexcept {
    return StrokeWidth_;
  }

  /// @brief Get the first Y-index of the interior cells. 
  size_t BeginCellY() const noexcept {
    return StrokeWidth_;
  }

  /// @brief Get the first Z-index of the interior cells. 
  size_t BeginCellZ() const noexcept {
    return StrokeWidth_;
  }

  /// @brief Get the X-index of the cell, \
  ///   following the last X-index of the the interior cells. 
  size_t EndCellX() const noexcept {
    return StrokeWidth_ + NumCellsX_;
  }

  /// @brief Get the Y-index of the cell, \
  ///   following the last Y-index of the the interior cells. 
  size_t EndCellY() const noexcept {
    return StrokeWidth_ + NumCellsY_;
  }

  /// @brief Get the Z-index of the cell, \
  ///    following the last Z-index of the the interior cells. 
  size_t EndCellZ() const noexcept {
    return StrokeWidth_ + NumCellsZ_;
  }

  /// @brief Get the volume of the cell at the 3D index.  
  real_t CellVolume(size_t ix, size_t iy, size_t iz) const noexcept {
    return CellVolume_(ix, iy, iz);
  }

  /// @brief Get the center position of the cell at the 3D index.  
  Vec3D<real_t> CellCenter(size_t ix, size_t iy, size_t iz) const noexcept {
    return CellCenter_(ix, iy, iz);
  }

  /// @brief Get the first X-index of the nodes. 
  size_t BeginNodeX() const noexcept {
    return StrokeWidth_;
  }

  /// @brief Get the first Y-index of the nodes. 
  size_t BeginNodeY() const noexcept {
    return StrokeWidth_;
  }

  /// @brief Get the first Z-index of the nodes. 
  size_t BeginNodeZ() const noexcept {
    return StrokeWidth_;
  }

  /// @brief Get the X-index of the node, \
  ///   following the last X-index of the the nodes. 
  size_t EndNodeX() const noexcept {
    return StrokeWidth_ + NumCellsX_ + 1;
  }

  /// @brief Get the Y-index of the node, \
  ///   following the last Y-index of the the nodes. 
  size_t EndNodeY() const noexcept {
    return StrokeWidth_ + NumCellsY_ + 1;
  }

  /// @brief Get the Z-index of the node, \
  ///   following the last Z-index of the the nodes. 
  size_t EndNodeZ() const noexcept {
    return StrokeWidth_ + NumCellsZ_ + 1;
  }

  /// @brief Get the position of the node at the 3D index.  
  Vec3D<real_t> NodePos(size_t ix, size_t iy, size_t iz) const noexcept {
    return NodePos_(ix, iy, iz);
  }

  /// @{
  real_t FaceAreaX(size_t ix, size_t iy, size_t iz) const noexcept {
    return FaceAreaX_(ix, iy, iz);
  }
  real_t FaceAreaY(size_t ix, size_t iy, size_t iz) const noexcept {
    return FaceAreaY_(ix, iy, iz);
  }
  real_t FaceAreaZ(size_t ix, size_t iy, size_t iz) const noexcept {
    return FaceAreaZ_(ix, iy, iz);
  }
  /// @}

  /// @{
  Vec3D<real_t> FaceNormalX(size_t ix, size_t iy, size_t iz) const noexcept {
    return FaceNormalX_(ix, iy, iz);
  }
  Vec3D<real_t> FaceNormalY(size_t ix, size_t iy, size_t iz) const noexcept {
    return FaceNormalY_(ix, iy, iz);
  }
  Vec3D<real_t> FaceNormalZ(size_t ix, size_t iy, size_t iz) const noexcept {
    return FaceNormalZ_(ix, iy, iz);
  }
  /// @}

}; // class Grid3D

class Quadrangle {
public:
  Quadrangle(std::initializer_list<Vec3D<real_t> const*>) {}
  real_t Area() const { return 0.0; }
  Vec3D<real_t> Normal() const { return {}; }
};

class Hexahedron {
public:
  Hexahedron(std::initializer_list<Vec3D<real_t> const*>) {}
  real_t Volume() const { return 0.0; }
  Vec3D<real_t> Center() const { return {}; }
};

void Grid3D::Build(size_t numGhostCells,
                   size_t numCellsX, size_t numCellsY, size_t numCellsZ,
                   std::function<Vec3D<real_t>(Vec3D<size_t>)> const& nodeTransform) {
  
  StormAssert(numGhostCells >= 1);
  StormAssert(numCellsX > 0 && numCellsY > 0 && numCellsZ > 0);

  // ----------------------
  // Assign the sizes.
  // ----------------------
  StrokeWidth_ = numGhostCells;
  NumCellsX_ = numCellsX, NumCellsY_ = numCellsY, NumCellsZ_ = numCellsZ;

  // ----------------------
  // Allocate the storage.
  // ----------------------
  size_t const sizeX = NumCellsX_ + 2*StrokeWidth_,
               sizeY = NumCellsY_ + 2*StrokeWidth_,
               sizeZ = NumCellsY_ + 2*StrokeWidth_;
  NodePos_.Assign(sizeX + 1, sizeY + 1, sizeZ + 1);
  CellVolume_.Assign(sizeX, sizeY, sizeZ);
  CellCenter_.Assign(sizeX, sizeY, sizeZ);
  FaceAreaX_.Assign(sizeX + 1, sizeY, sizeZ);
  FaceAreaY_.Assign(sizeX, sizeY + 1, sizeZ);
  FaceAreaZ_.Assign(sizeX, sizeY, sizeZ + 1);
  FaceNormalX_.Assign(sizeX + 1, sizeY, sizeZ);
  FaceNormalY_.Assign(sizeX, sizeY + 1, sizeZ);
  FaceNormalZ_.Assign(sizeX, sizeY, sizeZ + 1);

  // ----------------------
  // Build the node positions.
  // ----------------------
  for (size_t ix = 0; ix <= sizeX; ++ix) {
    for (size_t iy = 0; iy <= sizeY; ++iy) {
      for (size_t iz = 0; iz <= sizeZ; ++iz) {
        NodePos_(ix, iy, iz) = nodeTransform({ix, iy, iz});
      }
    }
  }

  // ----------------------
  // Compute the cell properties.
  // ----------------------
  for (size_t ix = 0; ix < sizeX; ++ix) {
    for (size_t iy = 0; iy < sizeY; ++iy) {
      for (size_t iz = 0; iz < sizeZ; ++iz) {
        Hexahedron const cell{
          &NodePos_(ix, iy, iz), &NodePos_(ix, iy, iz + 1),
          &NodePos_(ix, iy + 1, iz), &NodePos_(ix, iy + 1, iz + 1),
          &NodePos_(ix + 1, iy, iz), &NodePos_(ix + 1, iy, iz + 1),
          &NodePos_(ix + 1, iy + 1, iz), &NodePos_(ix + 1, iy + 1, iz + 1)
        };
        CellVolume_(ix, iy, iz) = 0.1*0.1*0.1; //cell.Volume();
        CellCenter_(ix, iy, iz) = {0.1*(ix + 0.5), 0.1*(iy + 0.5), 0.1*(iz + 0.5)}; //cell.Center();
      }
    }
  }

  // ----------------------
  // Compute the face properties.
  // ----------------------
  for (size_t ix = 0; ix <= sizeX; ++ix) {
    for (size_t iy = 0; iy < sizeY; ++iy) {
      for (size_t iz = 0; iz < sizeZ; ++iz) {
        Quadrangle const face{
          &NodePos_(ix, iy, iz), &NodePos_(ix, iy, iz + 1),
          &NodePos_(ix, iy + 1, iz), &NodePos_(ix, iy + 1, iz + 1)
        };
        FaceAreaX_(ix, iy, iz) = 0.1*0.1; //face.Area();
        FaceNormalX_(ix, iy, iz) = {1.0, 0.0, 0.0}; // face.Normal();
      }
    }
  }
  for (size_t ix = 0; ix < sizeX; ++ix) {
    for (size_t iy = 0; iy <= sizeY; ++iy) {
      for (size_t iz = 0; iz < sizeZ; ++iz) {
        Quadrangle const face{
          &NodePos_(ix, iy, iz), &NodePos_(ix, iy, iz + 1),
          &NodePos_(ix + 1, iy, iz), &NodePos_(ix + 1, iy, iz + 1)
        };
        FaceAreaY_(ix, iy, iz) = 0.1*0.1; //face.Area();
        FaceNormalY_(ix, iy, iz) = {0.0, 1.0, 0.0}; // face.Normal();
      }
    }
  }
  for (size_t ix = 0; ix < sizeX; ++ix) {
    for (size_t iy = 0; iy < sizeY; ++iy) {
      for (size_t iz = 0; iz <= sizeZ; ++iz) {
        Quadrangle const face{
          &NodePos_(ix, iy, iz), &NodePos_(ix, iy + 1, iz),
          &NodePos_(ix + 1, iy, iz), &NodePos_(ix + 1, iy + 1, iz)
        };
        FaceAreaZ_(ix, iy, iz) = 0.1*0.1; //face.Area();
        FaceNormalZ_(ix, iy, iz) = {0.0, 0.0, 1.0}; // face.Normal();
      }
    }
  }

} // Grid3D::Build

template<size_t Size>
void Print(std::ostream& out, 
           Grid3D const& G, 
           Array3D<Vec<real_t, Size>> const& Q) {

  out << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">" << std::endl;
  out << "<StructuredGrid WholeExtent=\"" << 
    G.BeginCellX() << " " << G.EndCellX() << " " <<
    G.BeginCellY() << " " << G.EndCellY() << " " <<
    G.BeginCellZ() << " " << G.EndCellZ() << "\">" << std::endl;
  out << "<Piece Extent=\"" << 
    G.BeginCellX() << " " << G.EndCellX() << " " <<
    G.BeginCellY() << " " << G.EndCellY() << " " <<
    G.BeginCellZ() << " " << G.EndCellZ() << "\">" << std::endl;

  // ----------------------
  // Print points.
  // ----------------------
  out << "<Points>" << std::endl;
  out << "<DataArray NumberOfComponents=\"3\" type=\"Float64\">" << std::endl;
  for (size_t iz = G.BeginNodeZ(); iz != G.EndNodeZ(); ++iz) {
    for (size_t iy = G.BeginNodeY(); iy != G.EndNodeY(); ++iy) {
      for (size_t ix = G.BeginNodeX(); ix != G.EndNodeX(); ++ix) {
        auto const nodePos = G.NodePos(ix, iy, iz);
        out << nodePos(0) << " " << nodePos(1) << " " << nodePos(2) << std::endl;
      }
    }
  }
  out << "</DataArray>" << std::endl;
  out << "</Points>" << std::endl;

  // ----------------------
  // Print arrays.
  // ----------------------
  out << "<CellData>" << std::endl;
  out << "<DataArray Name=\"rho\" NumberOfComponents=\"1\" type=\"Float64\">" << std::endl;
  for (size_t iz = G.BeginCellZ(); iz != G.EndCellZ(); ++iz) {
    for (size_t iy = G.BeginCellY(); iy != G.EndCellY(); ++iy) {
      for (size_t ix = G.BeginCellX(); ix != G.EndCellX(); ++ix) {
        out << Q(ix, iy, iz)(0) << std::endl;
      }
    }
  }
  out << "</DataArray>" << std::endl;
  out << "<DataArray Name=\"vel\" NumberOfComponents=\"3\" type=\"Float64\">" << std::endl;
  for (size_t iz = G.BeginCellZ(); iz != G.EndCellZ(); ++iz) {
    for (size_t iy = G.BeginCellY(); iy != G.EndCellY(); ++iy) {
      for (size_t ix = G.BeginCellX(); ix != G.EndCellX(); ++ix) {
        out << Q(ix, iy, iz)(2)/Q(ix, iy, iz)(0) << " " << 
               Q(ix, iy, iz)(3)/Q(ix, iy, iz)(0) << " " << 
               Q(ix, iy, iz)(4)/Q(ix, iy, iz)(0) << std::endl;
      }
    }
  }
  out << "</DataArray>" << std::endl;
  out << "</CellData>" << std::endl;

  out << "</Piece>" << std::endl;
  out << "</StructuredGrid>" << std::endl;
  out << "</VTKFile>" << std::endl;

} // Grid3D::Print

} // namespace Storm::Turbo
