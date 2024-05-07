//
// Created by 高建鑫 on 2024/4/21.
//

#ifndef KITAEV_KITAEV_TPS_SAMPLE_Z2_UPDATE_H
#define KITAEV_KITAEV_TPS_SAMPLE_Z2_UPDATE_H

#include "qlpeps/algorithm/vmc_update/wave_function_component.h"    //WaveFunctionComponent
#include "qlpeps/two_dim_tn/tensor_network_2d/tensor_network_2d.h"

namespace qlpeps {
    template<typename TenElemT, typename QNT>
    class KitaevTPSSampleZ2Flip : public WaveFunctionComponent<TenElemT, QNT> {
        using WaveFunctionComponentT = WaveFunctionComponent<TenElemT, QNT>;
    public:
        TensorNetwork2D<TenElemT, QNT> tn;

        KitaevTPSSampleZ2Flip(const size_t rows, const size_t cols) : WaveFunctionComponentT(rows, cols),
                                                                               tn(rows, cols) {}

        KitaevTPSSampleZ2Flip(const SplitIndexTPS<TenElemT, QNT> &sitps, const Configuration &config)
                : WaveFunctionComponentT(config), tn(config.rows(), config.cols()) {
            tn = TensorNetwork2D<TenElemT, QNT>(sitps, config);
            tn.GrowBMPSForRow(0, this->trun_para);
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
            tn.GrowBMPSForRow(0, this->trun_para);
            tn.GrowFullBTen(RIGHT, 0, 2, true);
            tn.InitBTen(LEFT, 0);
            this->amplitude = tn.Trace({0, 0}, HORIZONTAL);
        }

        void MonteCarloSweepUpdate(const SplitIndexTPS<TenElemT, QNT> &sitps,
                                   std::uniform_real_distribution<double> &u_double,
                                   std::vector<double> &accept_rates) {
            size_t flip_accept_num = 0;
            tn.GenerateBMPSApproach(UP, this->trun_para);
            for (size_t row = 0; row < tn.rows(); row++) {
                tn.InitBTen(LEFT, row);
                tn.GrowFullBTen(RIGHT, row, 2, true);
                for (size_t col = 0; col < tn.cols() - 1; col++) {
                    flip_accept_num += Z2FlipUpdate_({row, col}, {row, col + 1}, HORIZONTAL, sitps, u_double);
                    if (col < tn.cols() - 2) {
                        tn.BTenMoveStep(RIGHT);
                    }
                }
                if (row < tn.rows() - 1) {
                    tn.BMPSMoveStep(DOWN, this->trun_para);
                }
            }

            tn.DeleteInnerBMPS(LEFT);
            tn.DeleteInnerBMPS(RIGHT);

            tn.GenerateBMPSApproach(LEFT, this->trun_para);
            for (size_t col = 0; col < tn.cols(); col++) {
                tn.InitBTen(UP, col);
                tn.GrowFullBTen(DOWN, col, 2, true);
                for (size_t row = 0; row < tn.rows() - 1; row++) {
                    flip_accept_num += Z2FlipUpdate_({row, col}, {row + 1, col}, VERTICAL, sitps, u_double);
                    if (row < tn.rows() - 2) {
                        tn.BTenMoveStep(DOWN);
                    }
                }
                if (col < tn.cols() - 1) {
                    tn.BMPSMoveStep(RIGHT, this->trun_para);
                }
            }

            tn.DeleteInnerBMPS(UP);
            double bond_num = tn.cols() * (tn.rows() - 1) + tn.rows() * (tn.cols() - 1);
            accept_rates = {double(flip_accept_num) / bond_num};
        }

    private:
        bool Z2FlipUpdate_(const SiteIdx &site1, const SiteIdx &site2, BondOrientation bond_dir,
                           const SplitIndexTPS<TenElemT, QNT> &sitps,
                           std::uniform_real_distribution<double> &u_double) {

            //assert(sitps(site1)[this->config(site1)].GetIndexes() == sitps(site1)[this->config(site2)].GetIndexes());
            TenElemT psi_b = tn.ReplaceNNSiteTrace(site1, site2, bond_dir,
                                                   sitps(site1)[1 - this->config(site1)],
                                                   sitps(site2)[1 - this->config(site2)]);
            bool exchange;
            TenElemT &psi_a = this->amplitude;
            if (std::fabs(psi_b) >= std::fabs(psi_a)) {
                exchange = true;
            } else {
                double div = std::fabs(psi_b) / std::fabs(psi_a);
                double P = div * div;
                if (u_double(random_engine) < P) {
                    exchange = true;
                } else {
                    exchange = false;
                    return exchange;
                }
            }

            this->config(site1) = 1 - this->config(site1);
            this->config(site2) = 1 - this->config(site2);
            tn.UpdateSiteConfig(site1, this->config(site1), sitps);
            tn.UpdateSiteConfig(site2, this->config(site2), sitps);
            this->amplitude = psi_b;
            return exchange;
        }
    }; //KitaevTPSSampleZ2Flip

}//qlpeps


#endif //KITAEV_KITAEV_TPS_SAMPLE_Z2_UPDATE_H
