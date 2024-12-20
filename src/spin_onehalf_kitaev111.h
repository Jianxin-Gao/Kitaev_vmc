//
// Created by Jian-Xin Gao on 2024/9/20.
//

#ifndef KITAEV_SPIN_ONEHALF_KITAEV111_H
#define KITAEV_SPIN_ONEHALF_KITAEV111_H

#include "qlpeps/algorithm/vmc_update/model_energy_solver.h"      //ModelEnergySolver
#include "qlpeps/algorithm/vmc_update/model_measurement_solver.h" // ModelMeasurementSolver
#include "qlpeps/utility/helpers.h"
#include "qlcomplex.h"
#include "qlpeps/algorithm/vmc_update/wave_function_component.h"  //WaveFunctionComponent

namespace qlpeps {
using namespace qlten;

template<typename TenElemT, typename QNT>
class SpinOneHalfKitaev111 : public ModelEnergySolver<TenElemT, QNT>, ModelMeasurementSolver<TenElemT, QNT> {
  using SITPS = SplitIndexTPS<TenElemT, QNT>;
 public:
  SpinOneHalfKitaev111(void) = default;
  SpinOneHalfKitaev111(double Jx, double Jy, double Jz, double H) : Kx_(Jx), Ky_(Jy), Kz_(Jz), H_(H) {}

  using ModelEnergySolver<TenElemT, QNT>::ModelEnergySolver;

  template<typename WaveFunctionComponentType, bool calchols = true>
  TenElemT CalEnergyAndHoles(
      const SITPS *sitps,
      WaveFunctionComponentType *tps_sample,
      TensorNetwork2D<TenElemT, QNT> &hole_res);

  template<typename WaveFunctionComponentType>
  ObservablesLocal<TenElemT> SampleMeasure(
      const SITPS *sitps,
      WaveFunctionComponentType *tps_sample
  );
 private:
  double Kx_;
  double Ky_;
  double Kz_;
  double H_;
};

template<typename TenElemT, typename QNT>
template<typename WaveFunctionComponentType, bool calchols>
TenElemT SpinOneHalfKitaev111<TenElemT, QNT>::CalEnergyAndHoles(const SITPS *split_index_tps,
                                                                WaveFunctionComponentType *tps_sample,
                                                                TensorNetwork2D<TenElemT, QNT> &hole_res) {
  TenElemT energy(0);
  TensorNetwork2D<TenElemT, QNT> &tn = tps_sample->tn;
  const Configuration &config = tps_sample->config;
  const BMPSTruncatePara &trunc_para = WaveFunctionComponent<TenElemT, QNT>::trun_para.value();
  TenElemT inv_psi = 1.0 / (tps_sample->amplitude);
  std::vector<TenElemT> psi_gather;
  psi_gather.reserve(tn.rows() + tn.cols() - 2);
  tn.GenerateBMPSApproach(UP, trunc_para);

  for (size_t row = 0; row < tn.rows(); row++) {
    tn.InitBTen(LEFT, row);
    tn.GrowFullBTen(RIGHT, row, 1, true);
    // update the amplitude so that the error of ratio of amplitude can reduce by cancellation.
    tps_sample->amplitude = tn.Trace({row, 0}, HORIZONTAL);
    inv_psi = 1.0 / tps_sample->amplitude;
    psi_gather.push_back(tps_sample->amplitude);

    //Calculate horizontal bond energy contribution
    for (size_t col = 0; col < tn.cols(); col = col + 1) {
      const SiteIdx site1 = {row, col};
      //Calculate the holes
      if constexpr (calchols) {
        hole_res(site1) = Dag(tn.PunchHole(site1, HORIZONTAL));
      }
      size_t config1 = config(site1);  //  0 : up ; 1 : down
      energy += (config1 ? -H_ : H_); // H * sigma_z term
      TenElemT psi_flip = tn.ReplaceOneSiteTrace(site1, (*split_index_tps)(site1)[1 - config1], HORIZONTAL);
      TenElemT ratio = psi_flip * inv_psi;
      energy += ComplexConjugate(ratio) * H_ * (config1 ? std::complex<double>(1, -1) : std::complex<double>(1, 1));

      if (((row + col) & 1) == 0 && col + 1 < tn.cols()) {
        //Calculate horizontal bond energy contribution
        const SiteIdx site2 = {row, col + 1};

        TenElemT psi_ex = tn.ReplaceNNSiteTrace(site1, site2, HORIZONTAL,
                                                (*split_index_tps)(site1)[1 - config1],
                                                (*split_index_tps)(site2)[1 - config(site2)]);
        TenElemT ratio = psi_ex * inv_psi;
        energy += ComplexConjugate(ratio) * Kx_; // For Kx
      }

      if (col + 1 < tn.cols()) {
        tn.BTenMoveStep(RIGHT);
      }
    }

    if (row < tn.rows() - 1) {
      tn.BMPSMoveStep(DOWN, trunc_para);
    }
  }

  //Calculate vertical bond energy contribution
  tn.GenerateBMPSApproach(LEFT, trunc_para);
  for (size_t col = 0; col < tn.cols(); col++) {
    tn.InitBTen(UP, col);
    tn.GrowFullBTen(DOWN, col, 2, true);
    tps_sample->amplitude = tn.Trace({0, col}, VERTICAL);
    inv_psi = 1.0 / tps_sample->amplitude;
    psi_gather.push_back(tps_sample->amplitude);

    for (size_t row = 0; row + 1 < tn.rows(); row++) {
      const SiteIdx site1 = {row, col};
      const SiteIdx site2 = {row + 1, col};

      if ((row + col) % 2 == 0) {
        if (config(site1) == config(site2)) {
          energy += (Kz_); // For Kz
        } else {
          energy += -(Kz_); // For Kz
        }

      } else {
        TenElemT psi_ex = tn.ReplaceNNSiteTrace(site1, site2, VERTICAL,
                                                (*split_index_tps)(site1)[1 - config(site1)],
                                                (*split_index_tps)(site2)[1 - config(site2)]);
        TenElemT ratio = psi_ex * inv_psi;
        if (config(site1) == config(site2)) {
          energy += -ComplexConjugate(ratio) * (Ky_); // For Ky
        } else {
          energy += ComplexConjugate(ratio) * (Ky_); // For Ky
        }
      }
      if (row < tn.rows() - 2) {
        tn.BTenMoveStep(DOWN);
      }
    }

    if (col < tn.cols() - 1) {
      tn.BMPSMoveStep(RIGHT, trunc_para);
    }
  }
  return energy;
}

template<typename TenElemT, typename QNT>
template<typename WaveFunctionComponentType>
ObservablesLocal<TenElemT> SpinOneHalfKitaev111<TenElemT, QNT>::SampleMeasure(
    const SITPS *split_index_tps,
    WaveFunctionComponentType *tps_sample
) {
  ObservablesLocal<TenElemT> res;
  TenElemT energy(std::complex<double>(0.0, 0.0));
  TensorNetwork2D<TenElemT, QNT> &tn = tps_sample->tn;
  const size_t lx = tn.cols(), ly = tn.rows();
  res.bond_energys_loc.reserve(lx * ly * 2);
  res.two_point_functions_loc.reserve(lx / 2 * 3);
  const Configuration &config = tps_sample->config;
  const BMPSTruncatePara &trunc_para = WaveFunctionComponent<TenElemT, QNT>::trun_para.value();
  TenElemT inv_psi = 1.0 / (tps_sample->amplitude);
  std::vector<TenElemT> psi_gather;
  psi_gather.reserve(tn.rows() + tn.cols() - 2);
  tn.GenerateBMPSApproach(UP, trunc_para);
  res.one_point_functions_loc.reserve(
      tn.rows() * tn.cols() * 3); // 0:2:2*lx*ly-1 for <sigma_x>; 1:2:2*lx*ly-1 for <sigma_y>; 2*lx*ly:end for <sigma_z>
  for (size_t row = 0; row < ly; row++) {
    tn.InitBTen(LEFT, row);
    tn.GrowFullBTen(RIGHT, row, 1, true);
    // update the amplitude so that the error of ratio of amplitude can reduce by cancellation.
    tps_sample->amplitude = tn.Trace({row, 0}, HORIZONTAL);
    inv_psi = 1.0 / tps_sample->amplitude;
    psi_gather.push_back(tps_sample->amplitude);

    //Calculate horizontal bond energy contribution
    for (size_t col = 0; col < lx; col++) {
      std::complex<double> horizontal_bond_energy(0.0, 0.0);
      std::complex<double> Zeeman_energy(0.0, 0.0);
      const SiteIdx site1 = {row, col};
      size_t config1 = config(site1);  //  0 : up ; 1 : down
      Zeeman_energy += (config1 ? -H_ : H_); // H * sigma_z term
      TenElemT psi_flip = tn.ReplaceOneSiteTrace(site1, (*split_index_tps)(site1)[1 - config1], HORIZONTAL);
      TenElemT ratio = psi_flip * inv_psi;
      Zeeman_energy +=
          ComplexConjugate(ratio) * H_ * (config1 ? std::complex<double>(1, -1) : std::complex<double>(1, 1));
      res.one_point_functions_loc.push_back(ComplexConjugate(ratio)); // for <sigma_x>;
      res.one_point_functions_loc.push_back(ComplexConjugate(ratio)
                                                * (config1 ? std::complex<double>(0, -1) : std::complex<double>(0,
                                                                                                                1))); //for <sigma_y>

      energy += Zeeman_energy;

      if (((row + col) & 1) == 0 && col + 1 < tn.cols()) {
        //Calculate horizontal bond energy contribution
        const SiteIdx site2 = {row, col + 1};

        TenElemT psi_ex = tn.ReplaceNNSiteTrace(site1, site2, HORIZONTAL,
                                                (*split_index_tps)(site1)[1 - config1],
                                                (*split_index_tps)(site2)[1 - config(site2)]);
        ratio = psi_ex * inv_psi;
        horizontal_bond_energy += ComplexConjugate(ratio) * Kx_;
        energy += ComplexConjugate(ratio) * Kx_; // For Kx
        res.bond_energys_loc.push_back(horizontal_bond_energy);
      } else if (((row + col) & 1) != 0 && col + 1 < tn.cols()) {
        res.bond_energys_loc.push_back(horizontal_bond_energy);
      }

      if (col + 1 < tn.cols()) {
        tn.BTenMoveStep(RIGHT);
      }
    }

    if (row == tn.rows() / 2) { //measure correlation in the middle bonds
      SiteIdx site1 = {row, lx / 4};

      // sz(i) * sz(j)
      double sz1 = 1 - 2 * (double) config(site1);
      for (size_t i = 1; i <= lx / 2; i++) {
        SiteIdx site2 = {row, lx / 4 + i};
        double sz2 = 1 - 2 * (double) config(site2);
        res.two_point_functions_loc.push_back(sz1 * sz2);
      }

/*      std::vector<TenElemT> diag_corr(lx / 2);// sp(i) * sm(j) or sm(i) * sp(j), the valid channel
      tn(site1) = (*split_index_tps)(site1)[1 - config(site1)]; //temporally change
      tn.TruncateBTen(LEFT, lx / 4 + 1); // may be above two lines should be summarized as an API
      tn.GrowBTenStep(LEFT);//left boundary tensor just across Lx/4
      tn.GrowFullBTen(RIGHT, row, lx / 4 + 2, false); //environment for Lx/4 + 1 site
      for (size_t i = 1; i <= lx / 2; i++) {
        SiteIdx site2 = {row, lx / 4 + i};
        //sm(i) * sp(j)
        if (config(site2) == config(site1)) {
          diag_corr[i - 1] = 0.0;
        } else {
          TenElemT psi_ex = tn.ReplaceOneSiteTrace(site2, (*split_index_tps)(site2)[1 - config(site2)], HORIZONTAL);
          diag_corr[i - 1] = (psi_ex * inv_psi);
        }
        tn.BTenMoveStep(RIGHT);
      }
      tn(site1) = (*split_index_tps)(site1)[config(site1)]; // change back

      if (config(site1) == 1) {
        for (size_t i = 1; i <= lx / 2; i++) {  //sp(i) * sm(j) = 0
          res.two_point_functions_loc.push_back(0.0);
        }
        res.two_point_functions_loc.insert(res.two_point_functions_loc.end(), diag_corr.begin(), diag_corr.end());
      } else {
        res.two_point_functions_loc.insert(res.two_point_functions_loc.end(), diag_corr.begin(), diag_corr.end());
        for (size_t i = 1; i <= lx / 2; i++) {  //sm(i) * sp(j) = 0
          res.two_point_functions_loc.push_back(0.0);
        }
      }*/
    }
    if (row < tn.rows() - 1) {
      tn.BMPSMoveStep(DOWN, trunc_para);
    }
  }

  //Calculate vertical bond energy contribution
  tn.GenerateBMPSApproach(LEFT, trunc_para);
  for (size_t col = 0; col < tn.cols(); col++) {
    tn.InitBTen(UP, col);
    tn.GrowFullBTen(DOWN, col, 2, true);
    tps_sample->amplitude = tn.Trace({0, col}, VERTICAL);
    inv_psi = 1.0 / tps_sample->amplitude;
    psi_gather.push_back(tps_sample->amplitude);

    for (size_t row = 0; row < tn.rows() - 1; row++) {
      const SiteIdx site1 = {row, col};
      const SiteIdx site2 = {row + 1, col};
      std::complex<double> vertical_bond_energy(0.0, 0.0);
      if ((row + col) % 2 == 0) {
        if (config(site1) == config(site2)) {
          vertical_bond_energy += (Kz_);
          energy += (Kz_); // For Kz
        } else {
          vertical_bond_energy += -(Kz_);
          energy += -(Kz_); // For Kz
        }

      } else {
        TenElemT psi_ex = tn.ReplaceNNSiteTrace(site1, site2, VERTICAL,
                                                (*split_index_tps)(site1)[1 - config(site1)],
                                                (*split_index_tps)(site2)[1 - config(site2)]);
        TenElemT ratio = psi_ex * inv_psi;
        if (config(site1) == config(site2)) {
          vertical_bond_energy += -ComplexConjugate(ratio) * (Ky_);
          energy += -ComplexConjugate(ratio) * (Ky_); // For Ky
        } else {
          vertical_bond_energy += ComplexConjugate(ratio) * (Ky_);
          energy += ComplexConjugate(ratio) * (Ky_); // For Ky
        }
      }
      res.bond_energys_loc.push_back(vertical_bond_energy);
      if (row < tn.rows() - 2) {
        tn.BTenMoveStep(DOWN);
      }
    }
    if (col < tn.cols() - 1) {
      tn.BMPSMoveStep(RIGHT, trunc_para);
    }
  }
  res.energy_loc = energy;

  for (auto &spin_config : config) {
    res.one_point_functions_loc.push_back(1 - 2 * (double) spin_config); // For <sigma_z>
  }
  return res;
}

}//qlpeps






#endif //KITAEV_SPIN_ONEHALF_KITAEV111_H
