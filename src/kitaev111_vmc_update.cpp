
#include "./qlcomplex.h"
#include "qlpeps/algorithm/vmc_update/vmc_peps.h"
#include "qlpeps/algorithm/vmc_update/model_solvers/spin_onehalf_heisenberg_square.h"    // SpinOneHalfHeisenbergSquare
//#include "qlpeps/algorithm/vmc_update/wave_function_component_classes/square_tps_sample_3site_exchange.h"
#include "qlpeps/algorithm/vmc_update/wave_function_component_classes/square_tps_sample_full_space_nn_flip.h"
#include "params_parser.h"
#include "myutil.h"
#include "spin_onehalf_kitaev111.h"
//#include "Kitaev_tps_sample_Z2_update.h"
//#include "Kitaev_tps_sample_Z2_3site_update.h"

using namespace qlpeps;

using TPSSampleT = SquareTPSSampleFullSpaceNNFlip<TenElemT, U1QN>;

int main(int argc, char **argv) {
  boost::mpi::environment env;
  boost::mpi::communicator world;
  VMCUpdateParams params(argv[1]);

  qlten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);
  const size_t Ly = params.Ly;
  const size_t Lx = params.Lx;
  const size_t N = Lx * Ly;
  Configuration config_init(Ly, Lx);
  for (size_t row = 0; row < Ly; row++) {
    for (size_t col = 0; col < Lx; col++) {
      config_init({row, col}) = ((Ly) % 2 ? 0 : 1);
    }

  }
  qlpeps::VMCOptimizePara optimize_para(
      BMPSTruncatePara(params.Db_min, params.Db_max,
                       params.TruncErr,
                       params.MPSCompressScheme, 1e-8, 5),
      params.MC_samples, params.WarmUp,
      params.MCLocalUpdateSweepsBetweenSample,
      config_init,
      params.step_len,
      params.update_scheme,
      ConjugateGradientParams(params.CGMaxIter, params.CGTol, params.CGResidueRestart, params.CGDiagShift));

  using Model = SpinOneHalfKitaev111<TenElemT, U1QN>;
  VMCPEPSExecutor<TenElemT, U1QN, TPSSampleT, Model> *executor(nullptr);
  double Kx = params.Kx;
  double Ky = params.Ky;
  double Kz = params.Kz;
  Model kitaev111_solver(Kx, Ky, Kz, params.H);
  if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.qlten")) {// test if split index tps tensors exist
    executor = new VMCPEPSExecutor<TenElemT, U1QN, TPSSampleT, Model>(optimize_para,
                                                                      params.Ly, params.Lx,
                                                                      world, kitaev111_solver);
  } else {
    TPS<TenElemT, U1QN> tps = TPS<TenElemT, U1QN>(params.Ly, params.Lx);
    if (!tps.Load()) {
      std::cout << "Loading simple updated TPS files is broken." << std::endl;
      exit(-2);
    };
    executor = new VMCPEPSExecutor<TenElemT, U1QN, TPSSampleT, Model>(optimize_para, tps,
                                                                      world, kitaev111_solver);
  }
  executor->Execute();
  delete executor;

  return 0;
}
