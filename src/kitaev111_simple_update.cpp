/**
 * Created on 2024/06/05.
 *
 * Simple update for pure Kitaev model + [111] direction magnetic field. Complex number code.
 *
 * Hamiltonian:
 *     H = \sum_{<i,j>_\gamma} K_\gamma \sigma_i^\gamma \sigma_j^\gamma
 *          + h \sum_i (\sigma_i^x + \sigma_i^y +  \sigma_i^z)
 *
 * Same convention as in Ref https://www.nature.com/articles/s41467-022-28014-3
 */

#include "./kitaev_lattice_simple_update.h"
#include "./qlcomplex.h"
#include "./params_parser.h"

int main(int argc, char **argv) {
  SimpleUpdateParams params(argv[1]);
  // Get ham_hei_nn for x, y, z interaction
  Tensor ham_hei_nn_x = Tensor({pb_in, pb_out, pb_in, pb_out});
  Tensor ham_hei_nn_y = Tensor({pb_in, pb_out, pb_in, pb_out});
  Tensor ham_hei_nn_z = Tensor({pb_in, pb_out, pb_in, pb_out});

  auto Kx = params.Kx, Ky = params.Ky, Kz = params.Kz;
  auto H = params.H;
  // Kx * sigma_x * sigma_x + H * sigma_x * id + H * id * sigma_x;
  ham_hei_nn_x({0, 1, 1, 0}) = Kx;
  ham_hei_nn_x({1, 0, 1, 0}) = Kx;
  ham_hei_nn_x({0, 1, 0, 1}) = Kx;
  ham_hei_nn_x({1, 0, 0, 1}) = Kx;
  ham_hei_nn_x({0, 1, 0, 0}) = H; //sigma_x * id
  ham_hei_nn_x({1, 0, 0, 0}) = H;
  ham_hei_nn_x({0, 1, 1, 1}) = H;
  ham_hei_nn_x({1, 0, 1, 1}) = H;
  ham_hei_nn_x({0, 0, 0, 1}) = H; //id * sigma_x
  ham_hei_nn_x({0, 0, 1, 0}) = H;
  ham_hei_nn_x({1, 1, 0, 1}) = H;
  ham_hei_nn_x({1, 1, 1, 0}) = H;

  // Ky * sigma_y * simga_y + H * sigma_y * id + H * id * sigma_y;
  // sigma_y = ( 0, -i; i, 0)
  // 0 : up ; 1 : down
  ham_hei_nn_y({0, 1, 1, 0}) = Ky;
  ham_hei_nn_y({1, 0, 1, 0}) = -Ky;
  ham_hei_nn_y({0, 1, 0, 1}) = -Ky;
  ham_hei_nn_y({1, 0, 0, 1}) = Ky;
  ham_hei_nn_y({0, 1, 0, 0}) = TenElemT(0, H); //sigma_y * id, pure imaginary number
  ham_hei_nn_y({1, 0, 0, 0}) = TenElemT(0, -H);
  ham_hei_nn_y({0, 1, 1, 1}) = TenElemT(0, H);
  ham_hei_nn_y({1, 0, 1, 1}) = TenElemT(0, -H);
  ham_hei_nn_y({0, 0, 0, 1}) = TenElemT(0, H); //id * sigma_y
  ham_hei_nn_y({0, 0, 1, 0}) = TenElemT(0, -H);
  ham_hei_nn_y({1, 1, 0, 1}) = TenElemT(0, H);
  ham_hei_nn_y({1, 1, 1, 0}) = TenElemT(0, -H);

  // Kz * sigma_z * simga_z + H * sigma_z * id + H * id * sigma_z;
  ham_hei_nn_z({0, 0, 0, 0}) = Kz + 2 * H;
  ham_hei_nn_z({1, 1, 0, 0}) = -Kz;
  ham_hei_nn_z({0, 0, 1, 1}) = -Kz;
  ham_hei_nn_z({1, 1, 1, 1}) = Kz - 2 * H;

  qlten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);

  qlpeps::SimpleUpdatePara update_para(params.Step, params.Tau,
                                       params.Dmin, params.Dmax,
                                       params.TruncErr);

//  size_t Ly = params.Ly, Lx = params.Lx;
  qlpeps::SquareLatticePEPS<TenElemT, U1QN> peps0(pb_out, params.Ly, params.Lx);
  if (qlmps::IsPathExist(peps_path)) {
    peps0.Load(peps_path);
  } else {
    std::vector<std::vector<size_t>> activates(params.Ly, std::vector<size_t>(params.Lx));
    for (size_t y = 0; y < params.Ly; y++) {
      for (size_t x = 0; x < params.Lx; x++) {
        size_t sz_int = x + y;
        activates[y][x] = sz_int % 2;
      }
    }
    peps0.Initial(activates);
  }

  auto su_exe = new qlpeps::KitaevLatticeSimpleUpdateExecutor<TenElemT, U1QN>(update_para, peps0,
                                                                              ham_hei_nn_x,
                                                                              ham_hei_nn_y,
                                                                              ham_hei_nn_z);
  su_exe->Execute();
  auto tps = qlpeps::TPS<TenElemT, U1QN>(su_exe->GetPEPS());
  tps.Dump();
  su_exe->DumpResult(peps_path, true);

  return 0;
}