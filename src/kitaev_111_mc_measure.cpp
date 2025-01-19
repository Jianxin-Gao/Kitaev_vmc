//
// Created by 高建鑫 on 2024/9/20.
//

// SPDX-License-Identifier: LGPL-3.0-only

//
// Created by haoxinwang on 15/01/2024.
//


#include "qlpeps/algorithm/vmc_update/monte_carlo_measurement.h"
//#include "qlpeps/algorithm/vmc_update/vmc_optimize_para.h"
#include "spin_onehalf_kitaev111.h"
#include "qlpeps/algorithm/vmc_update/wave_function_component_classes/square_tps_sample_full_space_nn_flip.h"
//#include "./qldouble.h"
#include "qlcomplex.h"
#include "./params_parser.h"
#include "myutil.h"

using namespace qlpeps;

using TPSSampleT = SquareTPSSampleFullSpaceNNFlip<TenElemT, U1QN>;

int main(int argc, char **argv) {
    //boost::mpi::environment env;
    //boost::mpi::communicator world;
    MPI_Init(&argc, &argv);
    int rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    VMCUpdateParams params(argv[1]);

    qlten::hp_numeric::SetTensorManipulationThreads(params.ThreadNum);

    size_t N = params.Lx * params.Ly;
    qlpeps::MCMeasurementPara measurement_para(
            BMPSTruncatePara(params.Db_min, params.Db_max,
                             params.TruncErr,
                             params.MPSCompressScheme,
                             std::make_optional<double>(params.TruncErr),
                             std::make_optional<size_t>(10)),
            params.MC_samples, params.WarmUp,
            params.MCLocalUpdateSweepsBetweenSample,
            std::vector<size_t>{N / 2, N / 2},
            params.Ly, params.Lx);

    using Model = SpinOneHalfKitaev111<TenElemT, U1QN>;
    MonteCarloMeasurementExecutor<TenElemT, U1QN, TPSSampleT, Model> *executor(nullptr);
    double Kx = params.Kx;
    double Ky = params.Ky;
    double Kz = params.Kz;
    Model kitaev111_solver(Kx, Ky, Kz, params.H);
    if (IsFileExist(
            measurement_para.wavefunction_path + "/tps_ten0_0_0.qlten")) {// test if split index tps tensors exsit
        executor = new MonteCarloMeasurementExecutor<TenElemT, U1QN, TPSSampleT, Model>(measurement_para,
                                                                                        params.Ly, params.Lx,
                                                                                        MPI_COMM_WORLD,
                                                                                        kitaev111_solver);
    } else {
        TPS<QLTEN_Complex, U1QN> tps = TPS<QLTEN_Complex, U1QN>(params.Ly, params.Lx);
        if (!tps.Load()) {
            std::cout << "Loading simple updated TPS files is broken." << std::endl;
            exit(-2);
        };
        executor = new MonteCarloMeasurementExecutor<TenElemT, U1QN, TPSSampleT, Model>(measurement_para,
                                                                                        tps,
                                                                                        MPI_COMM_WORLD, kitaev111_solver);
    }

    executor->Execute();
    delete executor;
    //std::string bondinfo_filename = "energy_bonds" + std::to_string(params.Ly) + "-" + std::to_string(params.Lx);

    MPI_Finalize();
    return 0;
}
