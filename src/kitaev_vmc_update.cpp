//
// Created by 高建鑫 on 2024/4/14.
//

#include "./qldouble.h"
#include "qlpeps/algorithm/vmc_update/vmc_peps.h"
#include "qlpeps/algorithm/vmc_update/model_solvers/spin_onehalf_heisenberg_square.h"    // SpinOneHalfHeisenbergSquare
//#include "qlpeps/algorithm/vmc_update/wave_function_component_classes/square_tps_sample_3site_exchange.h"
//#include "qlpeps/algorithm/vmc_update/wave_function_component_classes/square_tps_sample_full_space_nn_flip.h"
#include "params_parser.h"
#include "myutil.h"
#include "spin_onehalf_kitaev.h"
//#include "Kitaev_tps_sample_Z2_update.h"
#include "kitaev_tps_sample_Z2_3site_update.h"

using namespace qlpeps;

using TPSSampleT = SquareTPSSample3SiteZ2Flip<TenElemT, U1QN>;

int main(int argc, char **argv) {
    boost::mpi::environment env;
    boost::mpi::communicator world;
    VMCUpdateParams params(argv[1]);

    qlten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);

    size_t N = params.Lx * params.Ly;
    qlpeps::VMCOptimizePara optimize_para(
            BMPSTruncatePara(params.Db_min, params.Db_max,
                             params.TruncErr,
                             params.MPSCompressScheme, 1e-8, 5),
            params.MC_samples, params.WarmUp,
            params.MCLocalUpdateSweepsBetweenSample,
            std::vector<size_t>{N / 2, N / 2},
            params.Ly, params.Lx,
            params.step_len,
            params.update_scheme,
            ConjugateGradientParams(params.CGMaxIter, params.CGTol, params.CGResidueRestart, params.CGDiagShift));


        using Model = SpinOneHalfKitaev<TenElemT, U1QN>;
        VMCPEPSExecutor<TenElemT, U1QN, TPSSampleT, Model> *executor(nullptr);
        double Jx = params.Kx;
        double Jy = params.Ky;
        double Jz = params.Kz;
        Model JxJyJzSolver(Jx, Jy, Jz);
        if (IsFileExist(optimize_para.wavefunction_path + "/tps_ten0_0_0.qlten")) {// test if split index tps tensors exist
            executor = new VMCPEPSExecutor<TenElemT, U1QN, TPSSampleT, Model>(optimize_para,
                                                                              params.Ly, params.Lx,
                                                                              world, JxJyJzSolver);
        } else {
            TPS<QLTEN_Double, U1QN> tps = TPS<QLTEN_Double, U1QN>(params.Ly, params.Lx);
            if (!tps.Load()) {
                std::cout << "Loading simple updated TPS files is broken." << std::endl;
                exit(-2);
            };
            executor = new VMCPEPSExecutor<TenElemT, U1QN, TPSSampleT, Model>(optimize_para, tps,
                                                                              world, JxJyJzSolver);
        }
        executor->Execute();
        delete executor;

    return 0;
}
