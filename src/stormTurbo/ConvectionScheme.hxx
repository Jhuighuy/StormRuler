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

#include <stormBase.hxx>
#include <stormSolvers/Mat.hxx>

#include <stormTurbo/Grid.hxx>
#include <stormTurbo/Array.hxx>
#include <stormTurbo/Parallel.hxx>

namespace Storm::Turbo {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract convection scheme.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value, size_t Size>
class ConvectionScheme3D {
public:

  virtual void Compute(Grid3D& G,
                       Array3D<Vec<Value, Size>>& Qhat,
                       Array3D<Vec<Value, Size>> const& Q) const = 0;

}; // class ConvectionScheme3D

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Upwind convection scheme.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value, size_t Size, 
         class FluxScheme>
class UpwindConvectionScheme3D final : 
  public ConvectionScheme3D<Value, Size> {
public:
  FluxScheme FluxScheme_;

public:

  void Compute(Grid3D& G,
               Array3D<Vec<Value, Size>>& Qhat,
               Array3D<Vec<Value, Size>> const& Q) const override;

}; // class UpwindConvectionScheme3D

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Piecewise-linear upwind convection scheme.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value, size_t Size, 
         class FluxScheme, class GradientScheme, class LimiterScheme>
class Upwind2ConvectionScheme3D final : 
  public ConvectionScheme3D<Value, Size> {
private:
  FluxScheme FluxScheme_;
  GradientScheme GradientScheme_;
  LimiterScheme LimiterScheme_;

public:

  void Compute(Grid3D& G,
               Array3D<Vec<Value, Size>>& Qhat,
               Array3D<Vec<Value, Size>> const& Q) const override;

}; // class Upwind2ConvectionScheme3D

template<class Value, size_t Size, 
         class FluxScheme>
void UpwindConvectionScheme3D<Value, Size, FluxScheme>::
                        Compute(Grid3D& G,
                                Array3D<Vec<Value, Size>>& Qhat,
                                Array3D<Vec<Value, Size>> const& Q) const {

  //Blas::Fill(Qhat, 0.0);

  // ----------------------
  // Evaluate the fluxes through the X-faces.
  // ----------------------
  ParFor(G.BeginCellY(), G.EndCellY(),
         G.BeginCellZ(), G.EndCellZ(), [&](size_t iy, size_t iz) {
    for (size_t ix = G.BeginCellX(); ix <= G.EndCellX(); ++ix) {
      Vec<Value, Size> const flux = 
        FluxScheme_(G.FaceNormalX(ix, iy, iz), 
                    Q(ix - 1, iy, iz), Q(ix, iy, iz));
      Qhat(ix, iy, iz) -= flux*G.FaceAreaX(ix, iy, iz);
      Qhat(ix - 1, iy, iz) += flux*G.FaceAreaX(ix, iy, iz);
    }
  });

  // ----------------------
  // Evaluate the fluxes through the Y-faces.
  // ----------------------
  ParFor(G.BeginCellX(), G.EndCellX(),
         G.BeginCellZ(), G.EndCellZ(), [&](size_t ix, size_t iz) {
    for (size_t iy = G.BeginCellY(); iy <= G.EndCellY(); ++iy) {
      Vec<Value, Size> const flux = 
        FluxScheme_(G.FaceNormalY(ix, iy, iz), 
                    Q(ix, iy, iz), Q(ix, iy - 1, iz));
      Qhat(ix, iy, iz) -= flux*G.FaceAreaY(ix, iy, iz);
      Qhat(ix, iy - 1, iz) += flux*G.FaceAreaY(ix, iy, iz);
    }
  });

  // ----------------------
  // Evaluate the fluxes through the Z-faces.
  // ----------------------
  ParFor(G.BeginCellX(), G.EndCellX(),
         G.BeginCellY(), G.EndCellY(), [&](size_t ix, size_t iy) {
    for (size_t iz = G.BeginCellZ(); iz <= G.EndCellZ(); ++iz) {
      Vec<Value, Size> const flux = 
        FluxScheme_(G.FaceNormalZ(ix, iy, iz), 
                    Q(ix, iy, iz), Q(ix, iy, iz - 1));
      Qhat(ix, iy, iz) -= flux*G.FaceAreaZ(ix, iy, iz);
      Qhat(ix, iy, iz - 1) += flux*G.FaceAreaZ(ix, iy, iz);
    }
  });

  // ----------------------
  // Finalize the fluxes summation.
  // ----------------------
  ParFor(G.BeginCellX(), G.EndCellX(),
         G.BeginCellY(), G.EndCellY(),
         G.BeginCellZ(), G.EndCellZ(), 
         [&] (size_t ix, size_t iy, size_t iz) {
   Qhat(ix, iy, iz) /= G.CellVolume(ix, iy, iz);
  });

} // UpwindConvectionScheme3D::Compute

} // namespace Storm::Turbo
