/*
 * File Name: params_case.h
 * Description: Declare CaseParams class used set parameters by users
 * Created by Hao-Xin on 2024/09/21.
 *
 */


#ifndef KITAEV_VMC_SRC_DMRG_PARAMS_CASE_H
#define KITAEV_VMC_SRC_DMRG_PARAMS_CASE_H

#include "qlmps/case_params_parser.h"
using qlmps::CaseParamsParserBasic;

struct CaseParams : public CaseParamsParserBasic {
  CaseParams(const char *pf) : CaseParamsParserBasic(pf) {
//    Geometry = ParseStr("Geometry");
    Ly = ParseInt("Ly");
    Lx = ParseInt("Lx");
    Kx = ParseDouble("Kx");
    Ky = ParseDouble("Ky");
    Kz = ParseDouble("Kz");
    H = ParseDouble("H");
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
    Perturbation = ParseDouble("Perturbation");
    wavelength = ParseInt("wavelength");
    noise = ParseDoubleVec("noise");
//    SymmetryMode = ParseInt("SymmetryMode");
  }

  std::string Geometry;
  size_t Ly;
  size_t Lx;
  double Kx;
  double Ky;
  double Kz;
  double H;
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

#endif //KITAEV_VMC_SRC_DMRG_PARAMS_CASE_H
