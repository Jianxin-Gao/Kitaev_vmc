//
// Created by 高建鑫 on 2024/4/13.
//

#ifndef KITAEV_SPIN_ONEHALF_KITAEV_H
#define KITAEV_SPIN_ONEHALF_KITAEV_H


#include "qlpeps/algorithm/vmc_update/model_energy_solver.h"      //ModelEnergySolver
#include "qlpeps/algorithm/vmc_update/model_measurement_solver.h" // ModelMeasurementSolver

namespace qlpeps {
    using namespace qlten;

    template<typename TenElemT, typename QNT>
    class SpinOneHalfKitaev : public ModelEnergySolver<TenElemT, QNT>, ModelMeasurementSolver<TenElemT, QNT> {
        using SITPS = SplitIndexTPS<TenElemT, QNT>;
    public:
        SpinOneHalfKitaev(void) = delete;
        SpinOneHalfKitaev(double Jx, double Jy, double Jz) : Jx_(Jx), Jy_(Jy), Jz_(Jz) {}


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
        double Jx_;
        double Jy_;
        double Jz_;
    };

    template<typename TenElemT, typename QNT>
    template<typename WaveFunctionComponentType, bool calchols>
    TenElemT SpinOneHalfKitaev<TenElemT, QNT>::CalEnergyAndHoles(const SITPS *split_index_tps,
                                                                           WaveFunctionComponentType *tps_sample,
                                                                           TensorNetwork2D<TenElemT, QNT> &hole_res) {
        const double bond_energy_extremly_large = 1.0e5;
        TenElemT energy(0);
        bool has_unreasonable_bond_energy(false);
        TensorNetwork2D<TenElemT, QNT> &tn = tps_sample->tn;
        const Configuration &config = tps_sample->config;
        const BMPSTruncatePara &trunc_para = SquareTPSSampleNNExchange<TenElemT, QNT>::trun_para;
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
            if ((row & 1) == 0) {
                for (size_t col = 0; col < tn.cols(); col = col + 1) {
                    const SiteIdx site1 = {row, col};
                    //Calculate the holes
                    if constexpr (calchols) {
                        hole_res(site1) = Dag(tn.PunchHole(site1, HORIZONTAL));
                    }
                    if ((col & 1) == 0) {
                        if (col + 1 < tn.cols()) {
                            //Calculate horizontal bond energy contribution
                            const SiteIdx site2 = {row, col + 1};

                            TenElemT psi_ex = tn.ReplaceNNSiteTrace(site1, site2, HORIZONTAL,
                                                                    (*split_index_tps)(site1)[1 - config(site1)],
                                                                    (*split_index_tps)(site2)[1 - config(site2)]);
                            TenElemT ratio = psi_ex * inv_psi;
                            energy += ratio * (-Jx_); // For Jx
                        }
                    }
                    if (col + 1 < tn.cols()) {
                        tn.BTenMoveStep(RIGHT);
                    }
                }
            } else {
                for (size_t col = 0; col < tn.cols(); col = col + 1){
                    const SiteIdx site1 = {row, col};
                    //Calculate the holes
                    if constexpr (calchols) {
                        hole_res(site1) = Dag(tn.PunchHole(site1, HORIZONTAL));
                    }
                    if ((col & 1) == 1) {
                        if (col + 1 < tn.cols()){
                            //Calculate horizontal bond energy contribution
                            const SiteIdx site2 = {row, col + 1};
                            TenElemT psi_ex = tn.ReplaceNNSiteTrace(site1, site2, HORIZONTAL,
                                                                    (*split_index_tps)(site1)[1 - config(site1)],
                                                                    (*split_index_tps)(site2)[1 - config(site2)]);
                            TenElemT ratio = psi_ex * inv_psi;
                            energy += ratio * (-Jx_); // For Jx
                        }
                    }
                    if (col + 1 < tn.cols()){
                        tn.BTenMoveStep(RIGHT);
                    }
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

            if ((col & 1) == 0) {
                for (size_t row = 0; row + 1 < tn.rows(); row ++) {
                    const SiteIdx site1 = {row, col};
                    const SiteIdx site2 = {row + 1, col};

                    if ((row & 1) == 0) {
                        if (config(site1) == config(site2)) {
                            energy += (-Jz_); // For Jz
                        } else {
                            energy += -(-Jz_); // For Jz
                        }

                    } else {
                        TenElemT psi_ex = tn.ReplaceNNSiteTrace(site1, site2, VERTICAL,
                                                                (*split_index_tps)(site1)[1 - config(site1)],
                                                                (*split_index_tps)(site2)[1 - config(site2)]);
                        TenElemT ratio = psi_ex * inv_psi;

                        if (config(site1) == config(site2)) {
                            energy += -ratio * (-Jy_); // For Jy
                        } else {
                            energy += ratio * (-Jy_); // For Jy
                        }
                    }
                    if (row < tn.rows() - 2) {
                        tn.BTenMoveStep(DOWN);
                    }
                }

            }
            else {
                for (size_t row = 0; row + 1 < tn.rows(); row ++) {
                    const SiteIdx site1 = {row, col};
                    const SiteIdx site2 = {row + 1, col};

                    if ((row & 1) == 1) {
                        if (config(site1) == config(site2)) {
                            energy += (-Jz_); // For Jz
                        } else {
                            energy += -(-Jz_); // For Jz
                        }
                    } else {
                        TenElemT psi_ex = tn.ReplaceNNSiteTrace(site1, site2, VERTICAL,
                                                                (*split_index_tps)(site1)[1 - config(site1)],
                                                                (*split_index_tps)(site2)[1 - config(site2)]);
                        TenElemT ratio = psi_ex * inv_psi;

                        if (config(site1) == config(site2)) {
                            energy += -ratio * (-Jy_); // For Jy
                        } else {
                            energy += ratio * (-Jy_); // For Jy
                        }
                    }
                    if (row < tn.rows() - 2) {
                        tn.BTenMoveStep(DOWN);
                    }
                }
            }


            if (col < tn.cols() - 1) {
                tn.BMPSMoveStep(RIGHT, trunc_para);
            }
        }
        if (has_unreasonable_bond_energy) {
            std::cout << "wave function amplitude estimation :" << std::endl;
            for (const auto &element : psi_gather) {
                std::cout << element << " ";
            }
            std::cout << std::endl;
        }
        return energy;
    }

    template<typename TenElemT, typename QNT>
    template<typename WaveFunctionComponentType>
    ObservablesLocal<TenElemT> SpinOneHalfKitaev<TenElemT, QNT>::SampleMeasure(
            const SITPS *split_index_tps,
            WaveFunctionComponentType *tps_sample
    ) {
        ObservablesLocal<TenElemT> res;
        TenElemT energy(0);
        const double bond_energy_extremly_large = 1.0e5;
        TensorNetwork2D<TenElemT, QNT> &tn = tps_sample->tn;
        const size_t lx = tn.cols(), ly = tn.rows();
        res.bond_energys_loc.reserve(lx * ly * 2);
        res.two_point_functions_loc.reserve(lx / 2 * 3);
        const Configuration &config = tps_sample->config;
        const BMPSTruncatePara &trunc_para = SquareTPSSampleNNExchange<TenElemT, QNT>::trun_para;
        TenElemT inv_psi = 1.0 / (tps_sample->amplitude);
        tn.GenerateBMPSApproach(UP, trunc_para);
        for (size_t row = 0; row < ly; row++) {
            tn.InitBTen(LEFT, row);
            tn.GrowFullBTen(RIGHT, row, 1, true);
            // update the amplitude so that the error of ratio of amplitude can reduce by cancellation.
            tps_sample->amplitude = tn.Trace({row, 0}, HORIZONTAL);
            inv_psi = 1.0 / tps_sample->amplitude;
            for (size_t col = 0; col < lx; col++) {
                const SiteIdx site1 = {row, col};
                if (col < tn.cols() - 1) {
                    //Calculate horizontal bond energy contribution
                    const SiteIdx site2 = {row, col + 1};
                    double horizontal_bond_energy;
                    if (config(site1) == config(site2)) {
                        horizontal_bond_energy = 0.25;
                    } else {
                        TenElemT psi_ex = tn.ReplaceNNSiteTrace(site1, site2, HORIZONTAL,
                                                                (*split_index_tps)(site1)[config(site2)],
                                                                (*split_index_tps)(site2)[config(site1)]);
                        TenElemT ratio = psi_ex * inv_psi;

                        horizontal_bond_energy = (-0.25 + ratio * 0.5);

                    }
                    energy += horizontal_bond_energy;
                    res.bond_energys_loc.push_back(horizontal_bond_energy);
                    tn.BTenMoveStep(RIGHT);
                }
            }
            if (row == tn.rows() / 2) { //measure correlation in the middle bonds
                SiteIdx site1 = {row, lx / 4};

                // sz(i) * sz(j)
                double sz1 = config(site1) - 0.5;
                for (size_t i = 1; i <= lx / 2; i++) {
                    SiteIdx site2 = {row, lx / 4 + i};
                    double sz2 = config(site2) - 0.5;
                    res.two_point_functions_loc.push_back(sz1 * sz2);
                }

                std::vector<TenElemT> diag_corr(lx / 2);// sp(i) * sm(j) or sm(i) * sp(j), the valid channel
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
            for (size_t row = 0; row < tn.rows() - 1; row++) {
                const SiteIdx site1 = {row, col};
                const SiteIdx site2 = {row + 1, col};
                double vertical_bond_energy;
                if (config(site1) == config(site2)) {
                    vertical_bond_energy = 0.25;
                } else {
                    TenElemT psi_ex = tn.ReplaceNNSiteTrace(site1, site2, VERTICAL,
                                                            (*split_index_tps)(site1)[config(site2)],
                                                            (*split_index_tps)(site2)[config(site1)]);
                    TenElemT ratio = psi_ex * inv_psi;

                    vertical_bond_energy = (-0.25 + ratio * 0.5);

                }
                energy += vertical_bond_energy;
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
        res.one_point_functions_loc.reserve(tn.rows() * tn.cols());
        for (auto &spin_config : config) {
            res.one_point_functions_loc.push_back((double) spin_config - 0.5);
        }
        return res;
    }

}//qlpeps






#endif //KITAEV_SPIN_ONEHALF_KITAEV_H
