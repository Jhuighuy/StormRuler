#include <iostream>
#include <fstream>

#include <stormTurbo/Grid.hxx>
#include <stormTurbo/Parallel.hxx>
#include <stormTurbo/Euler.hxx>
#include <stormTurbo/ConvectionScheme.hxx>

int main_turbo() {

  using namespace Storm;
  using namespace Storm::Turbo;

  std::cout << "Hello Turbo" << std::endl;

  Grid3D G;
  G.Build(1, 200, 50, 50, [](Vec3D<size_t> index) {
    //auto const z = index(0)/50.0;
    //auto const phi = std::atan(1.0)*index(1)/50.0; 
    //auto const r = 1.0 + index(2)/50.0;
    //return Vec3D<real_t>{ z, r*std::cos(phi), r*std::sin(phi) };
    auto const x = index(0)/50.0, y = index(1)/50.0, z = index(2)/50.0;
    return Vec3D<real_t>{ x, y, z };
  });

  Array3D<Vec<real_t, 5>> Q, F;

  Q.Assign(G.SizeX(), G.SizeY(), G.SizeZ());
  F.Assign(G.SizeX(), G.SizeY(), G.SizeZ());

  for (size_t ix = 0*G.BeginCellX(); ix != G.EndCellX(); ++ix) {
    for (size_t iy = 0*G.BeginCellY(); iy != G.EndCellY(); ++iy) {
      for (size_t iz = 0*G.BeginCellZ(); iz != G.EndCellZ(); ++iz) {
        Q(ix, iy, iz) = {10.0 + 0.0*std::sin(M_PI*0.1*ix), 1.0, 0.0, 0.0, 0.0};
        F(ix, iy, iz) = {0.0, 0.0, 0.0, 0.0, 0.0};
      }
    }
  }

  UpwindConvectionScheme3D<real_t, 5, decltype(&LaxFriedrichsFlux<real_t>)> scheme;
  scheme.FluxScheme_ = LaxFriedrichsFlux<real_t>;
  scheme.Compute(G, F, Q);

  std::ofstream out("out/grid.vts");
  Print(out, G, F);

  return 0;

}
