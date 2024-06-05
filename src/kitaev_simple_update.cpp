/**
 * Created on 2024/3/27.
 *
 * Simple update for pure Kitaev model. Real number code.
 *
 * Hamiltonian:
 *     H = \sum_{<i,j>_\gamma} K_\gamma \sigma_i^\gamma \sigma_j^\gamma
 *
 */

#include "./kitaev_lattice_simple_update.h"
#include "./qldouble.h"
#include "./params_parser.h"

int main(int argc, char **argv) {
  SimpleUpdateParams params(argv[1]);
  // Get ham_hei_nn for x, y, z interaction
  Tensor ham_hei_nn_x = Tensor({pb_in, pb_out, pb_in, pb_out});
  Tensor ham_hei_nn_y = Tensor({pb_in, pb_out, pb_in, pb_out});
  Tensor ham_hei_nn_z = Tensor({pb_in, pb_out, pb_in, pb_out});

  ham_hei_nn_x({0, 1, 1, 0}) = 1;
  ham_hei_nn_x({1, 0, 1, 0}) = 1;
  ham_hei_nn_x({0, 1, 0, 1}) = 1;
  ham_hei_nn_x({1, 0, 0, 1}) = 1;

  ham_hei_nn_y({0, 1, 1, 0}) = 1;
  ham_hei_nn_y({1, 0, 1, 0}) = -1;
  ham_hei_nn_y({0, 1, 0, 1}) = -1;
  ham_hei_nn_y({1, 0, 0, 1}) = 1;

  ham_hei_nn_z({0, 0, 0, 0}) = 1;
  ham_hei_nn_z({1, 1, 0, 0}) = -1;
  ham_hei_nn_z({0, 0, 1, 1}) = -1;
  ham_hei_nn_z({1, 1, 1, 1}) = 1;

  Tensor ham_x = params.Kx * ham_hei_nn_x;
  Tensor ham_y = params.Ky * ham_hei_nn_y;
  Tensor ham_z = params.Kz * ham_hei_nn_z;

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
                                                                              ham_x,
                                                                              ham_y,
                                                                              ham_z);
  su_exe->Execute();
  auto tps = qlpeps::TPS<TenElemT, U1QN>(su_exe->GetPEPS());
  tps.Dump();
  su_exe->DumpResult(peps_path, true);

  return 0;
}