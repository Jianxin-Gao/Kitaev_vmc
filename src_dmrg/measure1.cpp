#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <ctime>
#include "hilbert_space.h"
#include "operator.h"
#include "params_case.h"
#include "myutil.h"
#include "my_measure.h"

using namespace qlmps;
using namespace qlten;
using namespace std;

int main(int argc, char *argv[]) {
  namespace mpi = boost::mpi;
  mpi::environment env;
  mpi::communicator world;

  CaseParams params(argv[1]);
  size_t Lx = params.Lx;
  size_t N = Lx * params.Ly;
  cout << "The total number of sites: " << N << endl;
  clock_t startTime, endTime;
  startTime = clock();

  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);

  OperatorInitial();
  const SiteVec<TenElemT, TrivialRepQN> sites(N, pb_out);

  using FiniteMPST = qlmps::FiniteMPS<TenElemT, TrivialRepQN>;
  FiniteMPST mps(sites);

  Timer one_site_timer("measure  one site operators");
  MeasureOneSiteOp(mps, kMpsPath, {sigma_z, sigma_x, sigma_y}, {"sigma_z", "sigma_x", "sigma_y"});

  cout << "measured one point function.<====" << endl;
  one_site_timer.PrintElapsed();

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}