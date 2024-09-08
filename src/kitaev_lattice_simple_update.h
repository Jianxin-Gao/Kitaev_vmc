//
// Created by Jian-Xin Gao on 2024/3/16.
//

#ifndef KITAEV_KITAEV_LATTICE_SIMPLE_UPDATE_H
#define KITAEV_KITAEV_LATTICE_SIMPLE_UPDATE_H

#include "qlpeps/algorithm/simple_update/simple_update.h"
#include <complex>

namespace qlpeps {

    using namespace qlten;

    template<typename TenElemT, typename QNT>
    class KitaevLatticeSimpleUpdateExecutor : public SimpleUpdateExecutor<TenElemT, QNT> {
        using Tensor = QLTensor<TenElemT, QNT>;
        using PEPST = SquareLatticePEPS<TenElemT, QNT>;
    public:
        KitaevLatticeSimpleUpdateExecutor(const SimpleUpdatePara &update_para,
                                            const PEPST &peps_initial,
                                            const Tensor &ham_nn_x,
                                            const Tensor &ham_nn_y,
                                            const Tensor &ham_nn_z) :
                SimpleUpdateExecutor<TenElemT, QNT>(update_para, peps_initial), ham_nn_x_(ham_nn_x), ham_nn_y_(ham_nn_y), ham_nn_z_(ham_nn_z) {}

    private:
        void SetEvolveGate_(void) override {
            evolve_gate_nn_x_ = TaylorExpMatrix(this->update_para.tau, ham_nn_x_);
            evolve_gate_nn_y_ = TaylorExpMatrix(this->update_para.tau, ham_nn_y_);
            evolve_gate_nn_z_ = TaylorExpMatrix(this->update_para.tau, ham_nn_z_);
        }

        double SimpleUpdateSweep_(void) override;

        Tensor ham_nn_x_;
        Tensor ham_nn_y_;
        Tensor ham_nn_z_;
        Tensor evolve_gate_nn_x_;
        Tensor evolve_gate_nn_y_;
        Tensor evolve_gate_nn_z_;
    };


    template<typename TenElemT, typename QNT>
    double KitaevLatticeSimpleUpdateExecutor<TenElemT, QNT>::SimpleUpdateSweep_(void) {
        Timer simple_update_sweep_timer("simple_update_sweep");
        SimpleUpdateTruncatePara para(this->update_para.Dmin, this->update_para.Dmax, this->update_para.Trunc_err);
        double norm = 1.0;
        double e0 = 0.0;

#ifdef QLPEPS_TIMING_MODE
        Timer vertical_nn_projection_x_timer("vertical_nn_projection_x");
#endif
        for (size_t col = 0; col < this->lx_ - 1; col ++){
            if ((col & 1) == 0) {
                for (size_t row = 0; row < this->ly_ - 1; row = row + 2) {
                    ProjectionRes<std::complex<double>> result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                    norm = result.norm;

                    //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                    e0 += -std::log(norm) / this->update_para.tau;
                }
            }
            else {
                for (size_t row = 1; row < this->ly_; row = row + 2) {
                    ProjectionRes<std::complex<double>> result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                    norm = result.norm;
                    //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                    e0 += -std::log(norm) / this->update_para.tau;
                }
            }
        }

#ifdef QLPEPS_TIMING_MODE
        vertical_nn_projection_x_timer.PrintElapsed();
  Timer horizontal_nn_projection_y_timer("horizontal_nn_projection_y");
#endif
        for (size_t col = 0; col < this->lx_; col ++){
            if ((col & 1) == 0) {
                for (size_t row = 1; row < this->ly_ - 1; row = row + 2) {
                    ProjectionRes<std::complex<double>> result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_y_, {row, col}, VERTICAL, para);
                    norm = result.norm;
                    //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_y_, {row, col}, VERTICAL, para);
                    e0 += -std::log(norm) / this->update_para.tau;
                }
            }
            else {
                for (size_t row = 0; row < this->ly_ - 1; row = row + 2) {
                    ProjectionRes<std::complex<double>> result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_y_, {row, col}, VERTICAL, para);
                    norm = result.norm;
                    //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_y_, {row, col}, VERTICAL, para);
                    e0 += -std::log(norm) / this->update_para.tau;
                }
            }
        }
#ifdef QLPEPS_TIMING_MODE
        vertical_nn_projection_y_timer.PrintElapsed();
  Timer horizontal_nn_projection_z_timer("horizontal_nn_projection_z");
#endif
        for (size_t col = 0; col < this->lx_; col ++){
            if ((col & 1) == 0) {
                for (size_t row = 0; row < this->ly_ - 1; row = row + 2) {
                    ProjectionRes<std::complex<double>> result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_z_, {row, col}, VERTICAL, para);
                    norm = result.norm;
                    //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_z_, {row, col}, VERTICAL, para);
                    e0 += -std::log(norm) / this->update_para.tau;
                }
            }
            else {
                for (size_t row = 1; row < this->ly_ - 1; row = row + 2) {
                    ProjectionRes<std::complex<double>> result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_z_, {row, col}, VERTICAL, para);
                    norm = result.norm;
                    //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_z_, {row, col}, VERTICAL, para);
                    e0 += -std::log(norm) / this->update_para.tau;
                }
            }
        }
#ifdef QLPEPS_TIMING_MODE
        horizontal_nn_projection_z_timer.PrintElapsed();
#endif
        double sweep_time = simple_update_sweep_timer.Elapsed();
        auto [dmin, dmax] = this->peps_.GetMinMaxBondDim();
        std::cout << "Estimated E0 =" << std::setw(15) << std::setprecision(kEnergyOutputPrecision) << std::fixed
                  << std::right << e0
                  << " Dmin/Dmax = " << std::setw(2) << std::right << dmin << "/" << std::setw(2) << std::left << dmax
                  << " SweepTime = " << std::setw(8) << sweep_time
                  << std::endl;

        return norm;
    }
}

#endif //KITAEV_KITAEV_LATTICE_SIMPLE_UPDATE_H
