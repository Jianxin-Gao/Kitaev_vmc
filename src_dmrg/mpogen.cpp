#include "hilbert_space.h"
#include "operator.h"
#include "params_case.h"
#include "qlmps/qlmps.h"

using namespace std;
using namespace qlten;
using namespace qlmps;

int main(int argc, char *argv[]) {
  CaseParams params(argv[1]);
  const size_t Lx = params.Lx;
  const size_t Ly = params.Ly;
  const size_t N = Lx * params.Ly;
  cout << "The total number of sites: " << N << endl;
  double Kx = params.Kx;
  double Ky = params.Ky;
  double Kz = params.Kz;
  double H = params.H;
  cout << "Model parameter: Kx = " << Kx
       << ",\t Ky = " << Ky
       << ",\t Kz = " << Kz
       << ",\t H = " << params.H << endl;
  clock_t startTime, endTime;
  startTime = clock();

  if (!IsPathExist(qlmps::kMpoPath)) {
    CreatPath(qlmps::kMpoPath);
  }
  OperatorInitial();
  const SiteVec<TenElemT, TrivialRepQN> sites = SiteVec<TenElemT, TrivialRepQN>(N, pb_out);
  qlmps::MPOGenerator<TenElemT, TrivialRepQN> mpo_gen(sites, qn0);

  for (size_t i = 0; i < N; i++) {
    mpo_gen.AddTerm(H, sigma_x + sigma_y + sigma_z, i);
  }

  // x-x bond
  for (size_t i = 0; i < N; i++) {
    size_t x = i / Ly;
    size_t y = i % Ly;
    if (x < Lx - 1 && (x + y) % 2 == 0) {
      mpo_gen.AddTerm(Kx, sigma_x, i, sigma_x, i + Ly);
    }
  }
  // y-y & z-z bond
  for (size_t i = 0; i < N; i++) {
    size_t x = i / Ly;
    size_t y = i % Ly;
    if (y < Ly - 1) {
      if ((x + y) % 2 == 0) {
        mpo_gen.AddTerm(Kz, sigma_z, i, sigma_z, i + 1);
      } else {
        mpo_gen.AddTerm(Ky, sigma_y, i, sigma_y, i + 1);
      }
    }
  }
  
  auto mpo = mpo_gen.Gen();
  cout << "FiniteMPO generated." << endl;

  for (size_t i = 0; i < mpo.size(); i++) {
    std::string filename = qlmps::kMpoPath + "/" +
        qlmps::kMpoTenBaseName + std::to_string(i)
        + "." + qlmps::kQLTenFileSuffix;
    mpo.DumpTen(i, filename);
  }

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}
