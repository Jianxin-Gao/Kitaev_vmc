
/*
    measure2.cpp
    for measure spin, charge, on-site pair, single-particle correlation function.
    memory optimized and parallel version.
    usage:
        mpirun -n 2*Ly ./measure2
    Optional arguments:
      --start=
      --end=
    Which are set as start=Lx/4, end = 3*Lx/4+2 by default
*/

#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <ctime>
#include "hilbert_space.h"
#include "operator.h"
#include "params_case.h"
#include "myutil.h"
#include "my_measure.h"
#include <boost/serialization/complex.hpp>

#include "boost/mpi.hpp"

using std::cout;
using std::endl;
using std::vector;
using FiniteMPST = qlmps::FiniteMPS<TenElemT, TrivialRepQN>;
using qlmps::SiteVec;
using qlmps::MeasureTwoSiteOp;
using qlten::Timer;
using qlmps::MeasureGroupTask;
using qlmps::kMpsPath;

int main(int argc, char *argv[]) {
  namespace mpi = boost::mpi;
  mpi::environment env;
  mpi::communicator world;
  if (argc == 1) {
    std::cout << "Usage: \n mpirun -np 2*Ly ./measure2 <params file> --start=<start x coor> --end=<end x coor>\n";
    return 0;
  } else if (argc == 2) {
    std::cout << "No start and end parameters. set start as Lx/4 and end as 3*Lx/4" << std::endl;
    std::cout
        << "The complete usage is: \n mpirun -np 2*Ly ./measure2 <params file> --start=<start x coor> --end=<end x coor>\n";
  }

  clock_t startTime, endTime;
  startTime = clock();

  size_t beginx;
  size_t endx;
  bool start_argument_has = ParserMeasureSite(argc, argv, beginx, endx);

  CaseParams params(argv[1]);

  size_t Lx = params.Lx, Ly = params.Ly;
  size_t N = Lx * Ly;
  if (GetNumofMps() != N) {
    std::cout << "The number of mps files are inconsistent with mps size!" << std::endl;
    exit(1);
  }

  if (!start_argument_has) {
    beginx = Lx / 4;
    endx = beginx + Lx / 2 + 2;
  }

  OperatorInitial();

  const SiteVec<TenElemT, TrivialRepQN> sites = SiteVec<TenElemT, TrivialRepQN>(N, pb_out);
  FiniteMPST mps(sites);
  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);

  Timer two_site_timer("measure two site operators");
  std::vector<MeasureGroupTask> measure_tasks;
  measure_tasks.reserve(N);

  for (size_t y = 0; y < (Ly); ++y) {
    auto site1 = beginx * (Ly) + y;
    std::vector<size_t> site2;
    site2.reserve((Ly) * (endx - beginx));
    for (size_t j = (beginx + 2) * (Ly); j < endx * (Ly); j++) {
      site2.push_back(j);
    }
    measure_tasks.push_back(MeasureGroupTask(site1, site2));
  }

  Timer two_site_measure_timer("measure spin_ structure factors");
  mps.Load();
  MeasureTwoSiteOp(mps, sigma_z, sigma_z,
                   measure_tasks, "zzsf", world);
  if (world.rank() == 0) {
    std::cout << "measured sz sz correlation." << std::endl;
  }
  world.barrier();
  MeasureTwoSiteOp(mps, sigma_x, sigma_x, measure_tasks, "xxsf", world);
  MeasureTwoSiteOp(mps, sigma_y, sigma_y, measure_tasks, "yysf", world);
  two_site_measure_timer.PrintElapsed();

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

  return 0;
}


