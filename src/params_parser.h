#ifndef KITAEV_VMCPEPS_PARAMS_PARSER_H
#define KITAEV_VMCPEPS_PARAMS_PARSER_H

#include "qlmps/case_params_parser.h"
#include "qlpeps/algorithm/vmc_update/vmc_peps.h"

struct SimpleUpdateParams : public qlmps::CaseParamsParserBasic {
  SimpleUpdateParams(const char *f) : CaseParamsParserBasic(f) {
    Lx = ParseInt("Lx");
    Ly = ParseInt("Ly");
    RemoveCorner = ParseBool("RemoveCorner");
    H = ParseDouble("H");
    Kx = ParseDouble("Kx");
    Ky = ParseDouble("Ky");
    Kz = ParseDouble("Kz");
    TruncErr = ParseDouble("TruncErr");
    Dmin = ParseInt("Dmin");
    Dmax = ParseInt("Dmax");
    Tau = ParseDouble("Tau");
    Step = ParseInt("Step");
    ThreadNum = ParseInt("ThreadNum");
  }

  size_t Ly;
  size_t Lx;
  bool RemoveCorner;
  double H; // [111] direction magnetic field strength
  double Kx;
  double Ky;
  double Kz;
  double TruncErr;
  size_t Dmin;
  size_t Dmax;
  double Tau;
  size_t Step;
  size_t ThreadNum;

};

struct VMCUpdateParams : public qlmps::CaseParamsParserBasic {
  VMCUpdateParams(const char *f) : CaseParamsParserBasic(f) {
    Lx = ParseInt("Lx");
    Ly = ParseInt("Ly");
    H = ParseDouble("H");
    Kx = ParseDouble("Kx");
    Ky = ParseDouble("Ky");
    Kz = ParseDouble("Kz");
    Db_min = ParseInt("Dbmps_min");
    Db_max = ParseInt("Dbmps_max");
    TruncErr = ParseDouble("TruncErr");
    MC_samples = ParseInt("MC_samples");
    WarmUp = ParseInt("WarmUp");
    MCLocalUpdateSweepsBetweenSample = ParseInt("MCLocalUpdateSweepsBetweenSample");
    CGMaxIter = ParseInt("CGMaxIter");
    CGTol = ParseDouble("CGTol");
    CGResidueRestart = ParseInt("CGResidueRestart");
    CGDiagShift = ParseDouble("CGDiagShift");
    ReplicaTest = ParseBool("ReplicaTest");
    MPSCompressScheme = static_cast<qlpeps::CompressMPSScheme>(ParseInt("MPSCompressScheme"));
    RemoveCorner = ParseBool("RemoveCorner");
    size_t update_times = ParseInt("UpdateNum");
    step_len = std::vector<double>(update_times);
    if (update_times > 0) {
      step_len[0] = ParseDouble("StepLengthFirst");
      double step_len_change = ParseDouble("StepLengthDecrease");
      for (size_t i = 1; i < update_times; i++) {
        step_len[i] = step_len[0] - i * step_len_change;
      }
    }
    update_scheme = (qlpeps::WAVEFUNCTION_UPDATE_SCHEME) ParseInt("UpdateScheme");
    ThreadNum = ParseInt("ThreadNum");
  }

  size_t Ly;
  size_t Lx;
  double H; // [111] direction magnetic field strength
  size_t Db_min;
  size_t Db_max;
  double TruncErr;
  size_t MC_samples;
  size_t WarmUp;
  size_t MCLocalUpdateSweepsBetweenSample;
  size_t CGMaxIter;
  double CGTol;
  int CGResidueRestart;
  double CGDiagShift;
  bool ReplicaTest;
  qlpeps::CompressMPSScheme MPSCompressScheme;
  bool RemoveCorner;
  qlpeps::WAVEFUNCTION_UPDATE_SCHEME update_scheme;
  std::vector<double> step_len;
  size_t ThreadNum;
  double Kx;
  double Ky;
  double Kz;
};

struct DMRGCaseParams : public qlmps::CaseParamsParserBasic {
  DMRGCaseParams(const char *pf) : CaseParamsParserBasic(pf) {
//    Geometry = ParseStr("Geometry");
    Ly = ParseInt("Ly");
    Lx = ParseInt("Lx");
    RemoveCorner = ParseBool("RemoveCorner");
//    J1 = ParseDouble("J1");
    J2 = ParseDouble("J2");
//    J3 = ParseDouble("J3");
//    Dzz = ParseDouble("Dzz");
    Sweeps = ParseInt("Sweeps");
    Dmin = ParseInt("Dmin");
    Dmax = ParseInt("Dmax");
    CutOff = ParseDouble("CutOff");
    LanczErr = ParseDouble("LanczErr");
    MaxLanczIter = ParseInt("MaxLanczIter");
//    tau = ParseDouble("tau");
//    M = ParseInt("M");
//    ConvergeTolerance = ParseDouble("ConvergeTolerance");
//    TaylorOrder = ParseInt("TaylorOrder");
//    TaylorErr = ParseDouble("TaylorErr");
    Threads = ParseInt("Threads");
//    Perturbation=ParseDouble("Perturbation");
//    wavelength = ParseInt("wavelength");
    noise = ParseDoubleVec("noise");
//    SymmetryMode = ParseInt("SymmetryMode");
  }

  std::string Geometry;
  size_t Ly;
  size_t Lx;
  bool RemoveCorner;
  double J1;
  double J2;
  double J3;
  double Dzz;
  size_t Sweeps;
  size_t Dmin;
  size_t Dmax;
  double CutOff;
  double LanczErr;
  size_t MaxLanczIter;
  double tau;
  size_t M;
  double ConvergeTolerance;
  size_t TaylorOrder;
  double TaylorErr;
  size_t Threads;
  double Perturbation;
  size_t wavelength;
  std::vector<double> noise;
  size_t SymmetryMode;//useless upto now
};

#endif //KITAEV_VMCPEPS_PARAMS_PARSER_H
