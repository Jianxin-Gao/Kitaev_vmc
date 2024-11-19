/*
* Author: Hao-Xin Wang<wanghaoxin1996@gmail.com>
* Creation Date: 2024-05-07
*
* Description: QuantumLiquids/PEPS project.
*/

#ifndef QLPEPS_VMC_PEPS_SQUARE_TPS_SAMPLE_Z2_3SITE_UPDATE_H
#define QLPEPS_VMC_PEPS_SQUARE_TPS_SAMPLE_Z2_3SITE_UPDATE_H

#include "qlpeps/algorithm/vmc_update/wave_function_component.h"    // WaveFunctionComponent
#include "qlpeps/two_dim_tn/tensor_network_2d/tensor_network_2d.h"
#include "qlpeps/monte_carlo_tools/non_detailed_balance_mcmc.h"     // NonDBMCMCStateUpdate

namespace qlpeps {
template<typename TenElemT, typename QNT>
class SquareTPSSample3SiteZ2Flip : public WaveFunctionComponent<TenElemT, QNT> {
  using WaveFunctionComponentT = WaveFunctionComponent<TenElemT, QNT>;
 public:
  TensorNetwork2D<TenElemT, QNT> tn;

  SquareTPSSample3SiteZ2Flip(const size_t rows, const size_t cols) : WaveFunctionComponentT(rows, cols),
                                                                     tn(rows, cols) {}

  SquareTPSSample3SiteZ2Flip(const SplitIndexTPS<TenElemT, QNT> &sitps, const Configuration &config)
      : WaveFunctionComponentT(config), tn(config.rows(), config.cols()) {
    tn = TensorNetwork2D<TenElemT, QNT>(sitps, config);
    tn.GrowBMPSForRow(0, this->trun_para.value());
    tn.GrowFullBTen(RIGHT, 0, 2, true);
    tn.InitBTen(LEFT, 0);
    this->amplitude = tn.Trace({0, 0}, HORIZONTAL);
  }

  /**
   * @param sitps
   * @param occupancy_num
   */
  void RandomInit(const SplitIndexTPS<TenElemT, QNT> &sitps,
                  const std::vector<size_t> &occupancy_num) {
    this->config.Random(occupancy_num);
    tn = TensorNetwork2D<TenElemT, QNT>(sitps, this->config);
    tn.GrowBMPSForRow(0, this->trun_para.value());
    tn.GrowFullBTen(RIGHT, 0, 2, true);
    tn.InitBTen(LEFT, 0);
    this->amplitude = tn.Trace({0, 0}, HORIZONTAL);
  }

  void MonteCarloSweepUpdate(const SplitIndexTPS<TenElemT, QNT> &sitps,
                             std::uniform_real_distribution<double> &u_double,
                             std::vector<double> &accept_rates) {
    size_t flip_accept_num = 0;
    tn.GenerateBMPSApproach(UP, this->trun_para.value());
    for (size_t row = 0; row < tn.rows(); row++) {
      tn.InitBTen(LEFT, row);
      tn.GrowFullBTen(RIGHT, row, 3, true);
      for (size_t col = 0; col < tn.cols() - 2; col++) {
        flip_accept_num += Z23SiteUpdate_({row, col}, {row, col + 1}, {row, col + 2}, HORIZONTAL, sitps);
        if (col < tn.cols() - 3) {
          tn.BTenMoveStep(RIGHT);
        }
      }
      if (row < tn.rows() - 1) {
        tn.BMPSMoveStep(DOWN, this->trun_para.value());
      }
    }

    tn.DeleteInnerBMPS(LEFT);
    tn.DeleteInnerBMPS(RIGHT);

    tn.GenerateBMPSApproach(LEFT, this->trun_para.value());
    for (size_t col = 0; col < tn.cols(); col++) {
      tn.InitBTen(UP, col);
      tn.GrowFullBTen(DOWN, col, 3, true);
      for (size_t row = 0; row < tn.rows() - 2; row++) {
        flip_accept_num += Z23SiteUpdate_({row, col}, {row + 1, col}, {row + 2, col}, VERTICAL, sitps);
        if (row < tn.rows() - 3) {
          tn.BTenMoveStep(DOWN);
        }
      }
      if (col < tn.cols() - 1) {
        tn.BMPSMoveStep(RIGHT, this->trun_para.value());
      }
    }

    tn.DeleteInnerBMPS(UP);
    double total_flip_num = tn.cols() * (tn.rows() - 2) + tn.rows() * (tn.cols() - 2);
    accept_rates = {double(flip_accept_num) / total_flip_num};
  }

 private:
  ///< NB! physical dim == 2
  bool Z23SiteUpdate_(const SiteIdx &site1, const SiteIdx &site2, const SiteIdx &site3,
                      BondOrientation bond_dir,
                      const SplitIndexTPS<TenElemT, QNT> &sitps) {
    std::vector<std::vector<size_t>> local_config_set(4);
    std::vector<size_t> original_config = {this->config(site1), this->config(site2), this->config(site3)};
    local_config_set[0] = original_config;
    local_config_set[1] = {1 - this->config(site1), 1 - this->config(site2), this->config(site3)};
    local_config_set[2] = {1 - this->config(site1), this->config(site2), 1 - this->config(site3)};
    local_config_set[3] = {this->config(site1), 1 - this->config(site2), 1 - this->config(site3)};
    std::sort(local_config_set.begin(), local_config_set.end());
    std::vector<TenElemT> psi_set(4);
    double psi_abs_max = 0.0;
    size_t init_state;
    for (size_t i = 0; i < 4; i++) {
      if (local_config_set[i] == original_config) {
        psi_set[i] = this->amplitude;
        init_state = i;
      } else {
        psi_set[i] = tn.ReplaceTNNSiteTrace(site1, bond_dir,
                                            sitps(site1)[local_config_set[i][0]],
                                            sitps(site2)[local_config_set[i][1]],
                                            sitps(site3)[local_config_set[i][2]]);
      }
      psi_abs_max = std::max(psi_abs_max, psi_set[i]);
    }
    std::vector<double> weights(4);
    for (size_t i = 0; i < 4; i++) {
      weights[i] = std::norm(psi_set[i] / psi_abs_max);
    }

    size_t final_state = NonDBMCMCStateUpdate(init_state, weights, random_engine);
    if (final_state == init_state) {
      return false;
    } else {
      const std::vector<size_t> &new_local_config = local_config_set[final_state];
      this->config(site1) = new_local_config.at(0);
      this->config(site2) = new_local_config.at(1);
      this->config(site3) = new_local_config.at(2);
      this->amplitude = psi_set[final_state];
    }
    tn.UpdateSiteConfig(site1, this->config(site1), sitps);
    tn.UpdateSiteConfig(site2, this->config(site2), sitps);
    tn.UpdateSiteConfig(site3, this->config(site3), sitps);
    return true;
  }
};

}//qlpeps

#endif //QLPEPS_VMC_PEPS_SQUARE_TPS_SAMPLE_Z2_3SITE_UPDATE_H
