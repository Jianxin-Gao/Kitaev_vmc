//
// Created by Jian-Xin Gao on 2024/9/19.
//

#ifndef KITAEV_KITAEV_LATTICE_111_SIMPLE_UPDATE_H
#define KITAEV_KITAEV_LATTICE_111_SIMPLE_UPDATE_H

#include "qlpeps/algorithm/simple_update/simple_update.h"
#include <complex>

namespace qlpeps {

    using namespace qlten;

    template<typename TenElemT, typename QNT>
    class KitaevLattice111SimpleUpdateExecutor : public SimpleUpdateExecutor<TenElemT, QNT> {
        using Tensor = QLTensor<TenElemT, QNT>;
        using PEPST = SquareLatticePEPS<TenElemT, QNT>;
        using ComplexProjectionRes = ProjectionRes<TenElemT>;
    public:
        KitaevLattice111SimpleUpdateExecutor(const SimpleUpdatePara &update_para,
                                            const PEPST &peps_initial,
                                            const Tensor &ham_nn_x,
                                            const Tensor &ham_nn_y,
                                            const Tensor &ham_nn_z,
                                            const double H,
                                            const double Kz) :
                SimpleUpdateExecutor<TenElemT, QNT>(update_para, peps_initial), ham_nn_x_(ham_nn_x),
                ham_nn_y_(ham_nn_y),
                ham_nn_z_(ham_nn_z),
                H_(H),
                Kz_(Kz){}

    private:
        void SetEvolveGate_(void) override {
            Tensor ham_H_22_x_ = ham_nn_x_; // Kx * sigma_x * sigma_x + H/2 * [(sigma_x + sigma_y + sigma_z) * Id + Id * (sigma_x + sigma_y + sigma_z)];
            ham_H_22_x_({0, 1, 0, 0}) = std::complex<double>(H_, H_) / 2.0; //(sigma_x + sigma_y)* id
            ham_H_22_x_({1, 0, 0, 0}) = std::complex<double>(H_, -H_) / 2.0;
            ham_H_22_x_({0, 1, 1, 1}) = std::complex<double>(H_, H_) / 2.0;
            ham_H_22_x_({1, 0, 1, 1}) = std::complex<double>(H_, -H_) / 2.0;
            ham_H_22_x_({0, 0, 0, 1}) = std::complex<double>(H_, H_) / 2.0; //id * (sigma_x + sigma_y)*
            ham_H_22_x_({0, 0, 1, 0}) = std::complex<double>(H_, -H_) / 2.0;
            ham_H_22_x_({1, 1, 0, 1}) = std::complex<double>(H_, H_) / 2.0;
            ham_H_22_x_({1, 1, 1, 0}) = std::complex<double>(H_, -H_) / 2.0;
            ham_H_22_x_({0, 0, 0, 0}) = H_;
            ham_H_22_x_({1, 1, 1, 1}) = -H_;
            evolve_gate_H_22_x_ = TaylorExpMatrix(this->update_para.tau, ham_H_22_x_);

            Tensor ham_H_33_x_ = ham_nn_x_; // Kx * sigma_x * sigma_x + H/3 * [(sigma_x + sigma_y + sigma_z) * Id + Id * (sigma_x + sigma_y + sigma_z)];
            ham_H_33_x_({0, 1, 0, 0}) = std::complex<double>(H_, H_) / 3.0; //(sigma_x + sigma_y)* id
            ham_H_33_x_({1, 0, 0, 0}) = std::complex<double>(H_, -H_) / 3.0;
            ham_H_33_x_({0, 1, 1, 1}) = std::complex<double>(H_, H_) / 3.0;
            ham_H_33_x_({1, 0, 1, 1}) = std::complex<double>(H_, -H_) / 3.0;
            ham_H_33_x_({0, 0, 0, 1}) = std::complex<double>(H_, H_) / 3.0; //id * (sigma_x + sigma_y)*
            ham_H_33_x_({0, 0, 1, 0}) = std::complex<double>(H_, -H_) / 3.0;
            ham_H_33_x_({1, 1, 0, 1}) = std::complex<double>(H_, H_) / 3.0;
            ham_H_33_x_({1, 1, 1, 0}) = std::complex<double>(H_, -H_) / 3.0;
            ham_H_33_x_({0, 0, 0, 0}) = 2 * H_/3;
            ham_H_33_x_({1, 1, 1, 1}) = -2 * H_/3;
            evolve_gate_H_33_x_ = TaylorExpMatrix(this->update_para.tau, ham_H_33_x_);

            Tensor ham_H_23_y_ = ham_nn_y_; // Ky * sigma_y * sigma_y + H * [1/2 * (sigma_x + sigma_y + sigma_z)*Id + 1/3 * Id*(sigma_x + sigma_y + sigma_z)];
            ham_H_23_y_({0, 0, 0, 0}) = 5 * H_/6;
            ham_H_23_y_({0, 1, 0, 0}) = std::complex<double>(H_, H_) / 2.0;
            ham_H_23_y_({1, 0, 0, 0}) = std::complex<double>(H_, -H_) /2.0;
            ham_H_23_y_({1, 1, 0, 0}) = -H_/6;
            ham_H_23_y_({0, 0, 1, 0}) = std::complex<double>(H_, -H_) /3.0;
            ham_H_23_y_({1, 1, 1, 0}) = std::complex<double>(H_, -H_) /3.0;
            ham_H_23_y_({0, 0, 0, 1}) = std::complex<double>(H_, H_) /3.0;
            ham_H_23_y_({1, 1, 0, 1}) = std::complex<double>(H_, H_) /3.0;
            ham_H_23_y_({0, 0, 1, 1}) = H_ /6;
            ham_H_23_y_({0, 1, 1, 1}) = std::complex<double>(H_, H_) /2.0;
            ham_H_23_y_({1, 0, 1, 1}) = std::complex<double>(H_, -H_) /2.0;
            ham_H_23_y_({1, 1, 1, 1}) = -5 * H_/6;
            evolve_gate_H_23_y_ = TaylorExpMatrix(this->update_para.tau, ham_H_23_y_);

            Tensor ham_H_32_y_ = ham_nn_y_; // Ky * sigma_y * sigma_y + H * [1/3 * (sigma_x + sigma_y + sigma_z)*Id + 1/2 * Id*(sigma_x + sigma_y + sigma_z)];
            ham_H_32_y_({0, 0, 0, 0}) = 5 * H_/6;
            ham_H_32_y_({0, 1, 0, 0}) = std::complex<double>(H_, H_) / 3.0;
            ham_H_32_y_({1, 0, 0, 0}) = std::complex<double>(H_, -H_) /3.0;
            ham_H_32_y_({1, 1, 0, 0}) = H_/6;
            ham_H_32_y_({0, 0, 1, 0}) = std::complex<double>(H_, -H_) /2.0;
            ham_H_32_y_({1, 1, 1, 0}) = std::complex<double>(H_, -H_) /2.0;
            ham_H_32_y_({0, 0, 0, 1}) = std::complex<double>(H_, H_) /2.0;
            ham_H_32_y_({1, 1, 0, 1}) = std::complex<double>(H_, H_) /2.0;
            ham_H_32_y_({0, 0, 1, 1}) = -H_ /6;
            ham_H_32_y_({0, 1, 1, 1}) = std::complex<double>(H_, H_) /3.0;
            ham_H_32_y_({1, 0, 1, 1}) = std::complex<double>(H_, -H_) /3.0;
            ham_H_32_y_({1, 1, 1, 1}) = -5 * H_/6;
            evolve_gate_H_32_y_ = TaylorExpMatrix(this->update_para.tau, ham_H_32_y_);

            Tensor ham_H_33_y_ = ham_nn_y_; // Ky * sigma_y * sigma_y + H * [1/3 * (sigma_x + sigma_y + sigma_z)*Id + 1/3 * Id*(sigma_x + sigma_y + sigma_z)];
            ham_H_33_y_({0, 1, 0, 0}) = std::complex<double>(H_, H_) / 3.0; //(sigma_x + sigma_y)* id
            ham_H_33_y_({1, 0, 0, 0}) = std::complex<double>(H_, -H_) / 3.0;
            ham_H_33_y_({0, 1, 1, 1}) = std::complex<double>(H_, H_) / 3.0;
            ham_H_33_y_({1, 0, 1, 1}) = std::complex<double>(H_, -H_) / 3.0;
            ham_H_33_y_({0, 0, 0, 1}) = std::complex<double>(H_, H_) / 3.0; //id * (sigma_x + sigma_y)*
            ham_H_33_y_({0, 0, 1, 0}) = std::complex<double>(H_, -H_) / 3.0;
            ham_H_33_y_({1, 1, 0, 1}) = std::complex<double>(H_, H_) / 3.0;
            ham_H_33_y_({1, 1, 1, 0}) = std::complex<double>(H_, -H_) / 3.0;
            ham_H_33_y_({0, 0, 0, 0}) = 2 * H_ / 3.0;  //( H * sigma_z * id + H * id * sigma_z)/3
            ham_H_33_y_({1, 1, 1, 1}) = -2 * H_ / 3.0;
            evolve_gate_H_33_y_ = TaylorExpMatrix(this->update_para.tau, ham_H_33_y_);

            Tensor ham_H_22_y_ = ham_nn_y_; // Ky * sigma_y * sigma_y + H * [1/2 * (sigma_x + sigma_y + sigma_z)*Id + 1/2 * Id*(sigma_x + sigma_y + sigma_z)];
            ham_H_22_y_({0, 1, 0, 0}) = std::complex<double>(H_, H_) / 2.0; //(sigma_x + sigma_y)* id
            ham_H_22_y_({1, 0, 0, 0}) = std::complex<double>(H_, -H_) / 2.0;
            ham_H_22_y_({0, 1, 1, 1}) = std::complex<double>(H_, H_) / 2.0;
            ham_H_22_y_({1, 0, 1, 1}) = std::complex<double>(H_, -H_) / 2.0;
            ham_H_22_y_({0, 0, 0, 1}) = std::complex<double>(H_, H_) / 2.0; //id * (sigma_x + sigma_y)*
            ham_H_22_y_({0, 0, 1, 0}) = std::complex<double>(H_, -H_) / 2.0;
            ham_H_22_y_({1, 1, 0, 1}) = std::complex<double>(H_, H_) / 2.0;
            ham_H_22_y_({1, 1, 1, 0}) = std::complex<double>(H_, -H_) / 2.0;
            ham_H_22_y_({0, 0, 0, 0}) =  H_ ;  //( H * sigma_z * id + H * id * sigma_z)/2
            ham_H_22_y_({1, 1, 1, 1}) = -H_ ;
            evolve_gate_H_22_y_ = TaylorExpMatrix(this->update_para.tau, ham_H_22_y_);

            Tensor ham_H_31_y_ = ham_nn_y_; // Ky * sigma_y * sigma_y + H * [1/3 * (sigma_x + sigma_y + sigma_z)*Id +  Id*(sigma_x + sigma_y + sigma_z)];
            ham_H_31_y_({0, 0, 0, 0}) = 4 * H_/3;
            ham_H_31_y_({0, 1, 0, 0}) = std::complex<double>(H_, H_) / 3.0;
            ham_H_31_y_({1, 0, 0, 0}) = std::complex<double>(H_, -H_) /3.0;
            ham_H_31_y_({1, 1, 0, 0}) = 2 * H_/3;
            ham_H_31_y_({0, 0, 1, 0}) = std::complex<double>(H_, -H_);
            ham_H_31_y_({1, 1, 1, 0}) = std::complex<double>(H_, -H_);
            ham_H_31_y_({0, 0, 0, 1}) = std::complex<double>(H_, H_);
            ham_H_31_y_({1, 1, 0, 1}) = std::complex<double>(H_, H_);
            ham_H_31_y_({0, 0, 1, 1}) = -2 * H_ /3;
            ham_H_31_y_({0, 1, 1, 1}) = std::complex<double>(H_, H_) /3.0;
            ham_H_31_y_({1, 0, 1, 1}) = std::complex<double>(H_, -H_) /3.0;
            ham_H_31_y_({1, 1, 1, 1}) = -4 * H_/3;
            evolve_gate_H_31_y_ = TaylorExpMatrix(this->update_para.tau, ham_H_31_y_);

            Tensor ham_H_22_z_ = ham_nn_z_; // Kz * sigma_z * sigma_z + H * [1/2 * (sigma_x + sigma_y + sigma_z)*Id +  1/2 * Id*(sigma_x + sigma_y + sigma_z)];
            ham_H_22_z_({0, 0, 0, 0}) = Kz_ +  H_ ;
            ham_H_22_z_({1, 1, 1, 1}) = Kz_ -  H_ ;
            ham_H_22_z_({0, 1, 0, 0}) = std::complex<double>(H_, H_) / 2.0; //(sigma_x + sigma_y)* id
            ham_H_22_z_({1, 0, 0, 0}) = std::complex<double>(H_, -H_) / 2.0;
            ham_H_22_z_({0, 1, 1, 1}) = std::complex<double>(H_, H_) / 2.0;
            ham_H_22_z_({1, 0, 1, 1}) = std::complex<double>(H_, -H_) /2.0;
            ham_H_22_z_({0, 0, 0, 1}) = std::complex<double>(H_, H_) / 2.0; //id * (sigma_x + sigma_y)*
            ham_H_22_z_({0, 0, 1, 0}) = std::complex<double>(H_, -H_) / 2.0;
            ham_H_22_z_({1, 1, 0, 1}) = std::complex<double>(H_, H_) / 2.0;
            ham_H_22_z_({1, 1, 1, 0}) = std::complex<double>(H_, -H_) / 2.0;
            evolve_gate_H_22_z_ = TaylorExpMatrix(this->update_para.tau, ham_H_22_z_);

            Tensor ham_H_33_z_ = ham_nn_z_; // Kz * sigma_z * sigma_z + H * [1/3 * (sigma_x + sigma_y + sigma_z)*Id +  1/3 * Id*(sigma_x + sigma_y + sigma_z)];
            ham_H_33_z_({0, 0, 0, 0}) = Kz_ +  2 * H_/3 ;
            ham_H_33_z_({1, 1, 1, 1}) = Kz_ -  2 * H_/3 ;
            ham_H_33_z_({0, 1, 0, 0}) = std::complex<double>(H_, H_) / 3.0; //(sigma_x + sigma_y)* id
            ham_H_33_z_({1, 0, 0, 0}) = std::complex<double>(H_, -H_) / 3.0;
            ham_H_33_z_({0, 1, 1, 1}) = std::complex<double>(H_, H_) / 3.0;
            ham_H_33_z_({1, 0, 1, 1}) = std::complex<double>(H_, -H_) /3.0;
            ham_H_33_z_({0, 0, 0, 1}) = std::complex<double>(H_, H_) / 3.0; //id * (sigma_x + sigma_y)*
            ham_H_33_z_({0, 0, 1, 0}) = std::complex<double>(H_, -H_) / 3.0;
            ham_H_33_z_({1, 1, 0, 1}) = std::complex<double>(H_, H_) / 3.0;
            ham_H_33_z_({1, 1, 1, 0}) = std::complex<double>(H_, -H_) / 3.0;
            evolve_gate_H_33_z_ = TaylorExpMatrix(this->update_para.tau, ham_H_33_z_);

            Tensor ham_H_32_z_ = ham_nn_z_; // Kz * sigma_z * sigma_z + H * [1/3 * (sigma_x + sigma_y + sigma_z)*Id +  1/2 * Id*(sigma_x + sigma_y + sigma_z)];
            ham_H_32_z_({0, 0, 0, 0}) = Kz_ +  5 * H_/6;
            ham_H_32_z_({0, 1, 0, 0}) = std::complex<double>(H_, H_) /3.0;
            ham_H_32_z_({1, 0, 0, 0}) = std::complex<double>(H_, -H_) /3.0 ;
            ham_H_32_z_({1, 1, 0, 0}) = -Kz_ + H_/6 ;
            ham_H_32_z_({0, 0, 1, 0}) = std::complex<double>(H_, -H_) /2.0;
            ham_H_32_z_({1, 1, 1, 0}) = std::complex<double>(H_, -H_) /2.0;
            ham_H_32_z_({0, 0, 0, 1}) = std::complex<double>(H_, H_) /2.0;
            ham_H_32_z_({1, 1, 0, 1}) = std::complex<double>(H_, H_) /2.0;
            ham_H_32_z_({0, 0, 1, 1}) = -Kz_ - H_/6;
            ham_H_32_z_({0, 1, 1, 1}) = std::complex<double>(H_, H_) /3.0;
            ham_H_32_z_({1, 0, 1, 1}) = std::complex<double>(H_, -H_) /3.0;
            ham_H_32_z_({1, 1, 1, 1}) = Kz_ - 5*H_/6 ;
            evolve_gate_H_32_z_ = TaylorExpMatrix(this->update_para.tau, ham_H_32_z_);

            Tensor ham_H_23_z_ = ham_nn_z_; // Kz * sigma_z * sigma_z + H * [1/2 * (sigma_x + sigma_y + sigma_z)*Id +  1/3 * Id*(sigma_x + sigma_y + sigma_z)];
            ham_H_23_z_({0, 0, 0, 0}) = Kz_ +  5 * H_/6;
            ham_H_23_z_({0, 1, 0, 0}) = std::complex<double>(H_, H_) /2.0;
            ham_H_23_z_({1, 0, 0, 0}) = std::complex<double>(H_, -H_) /2.0 ;
            ham_H_23_z_({1, 1, 0, 0}) = -Kz_ - H_/6 ;
            ham_H_23_z_({0, 0, 1, 0}) = std::complex<double>(H_, -H_) /3.0;
            ham_H_23_z_({1, 1, 1, 0}) = std::complex<double>(H_, -H_) /3.0;
            ham_H_23_z_({0, 0, 0, 1}) = std::complex<double>(H_, H_) /3.0;
            ham_H_23_z_({1, 1, 0, 1}) = std::complex<double>(H_, H_) /3.0;
            ham_H_23_z_({0, 0, 1, 1}) = -Kz_ + H_/6;
            ham_H_23_z_({0, 1, 1, 1}) = std::complex<double>(H_, H_) /2.0;
            ham_H_23_z_({1, 0, 1, 1}) = std::complex<double>(H_, -H_) /2.0;
            ham_H_23_z_({1, 1, 1, 1}) = Kz_ - 5*H_/6 ;
            evolve_gate_H_23_z_ = TaylorExpMatrix(this->update_para.tau, ham_H_23_z_);

            Tensor ham_H_31_z_ = ham_nn_z_; // Kz * sigma_z * sigma_z + H * [1/3 * (sigma_x + sigma_y + sigma_z)*Id +  Id*(sigma_x + sigma_y + sigma_z)];
            ham_H_31_z_({0, 0, 0, 0}) = Kz_ +  4 * H_/3;
            ham_H_31_z_({0, 1, 0, 0}) = std::complex<double>(H_, H_) /3.0;
            ham_H_31_z_({1, 0, 0, 0}) = std::complex<double>(H_, -H_) /3.0 ;
            ham_H_31_z_({1, 1, 0, 0}) = -Kz_ + 2*H_/3 ;
            ham_H_31_z_({0, 0, 1, 0}) = std::complex<double>(H_, -H_);
            ham_H_31_z_({1, 1, 1, 0}) = std::complex<double>(H_, -H_);
            ham_H_31_z_({0, 0, 0, 1}) = std::complex<double>(H_, H_);
            ham_H_31_z_({1, 1, 0, 1}) = std::complex<double>(H_, H_);
            ham_H_31_z_({0, 0, 1, 1}) = -Kz_ -2* H_/3;
            ham_H_31_z_({0, 1, 1, 1}) = std::complex<double>(H_, H_) /3.0;
            ham_H_31_z_({1, 0, 1, 1}) = std::complex<double>(H_, -H_) /3.0;
            ham_H_31_z_({1, 1, 1, 1}) = Kz_ - 4*H_/3 ;
            evolve_gate_H_31_z_ = TaylorExpMatrix(this->update_para.tau, ham_H_31_z_);

            Tensor ham_H_13_z_ = ham_nn_z_; // Kz * sigma_z * sigma_z + H * [ (sigma_x + sigma_y + sigma_z)*Id +  1/3 * Id*(sigma_x + sigma_y + sigma_z)];
            ham_H_13_z_({0, 0, 0, 0}) = Kz_ +  4 * H_/3;
            ham_H_13_z_({0, 1, 0, 0}) = std::complex<double>(H_, H_);
            ham_H_13_z_({1, 0, 0, 0}) = std::complex<double>(H_, -H_) ;
            ham_H_13_z_({1, 1, 0, 0}) = -Kz_ - 2*H_/3 ;
            ham_H_13_z_({0, 0, 1, 0}) = std::complex<double>(H_, -H_) /3.0;
            ham_H_13_z_({1, 1, 1, 0}) = std::complex<double>(H_, -H_) /3.0;
            ham_H_13_z_({0, 0, 0, 1}) = std::complex<double>(H_, H_) /3.0;
            ham_H_13_z_({1, 1, 0, 1}) = std::complex<double>(H_, H_) /3.0;
            ham_H_13_z_({0, 0, 1, 1}) = -Kz_ +2* H_/3;
            ham_H_13_z_({0, 1, 1, 1}) = std::complex<double>(H_, H_);
            ham_H_13_z_({1, 0, 1, 1}) = std::complex<double>(H_, -H_);
            ham_H_13_z_({1, 1, 1, 1}) = Kz_ - 4*H_/3 ;
            evolve_gate_H_13_z_ = TaylorExpMatrix(this->update_para.tau, ham_H_13_z_);

            evolve_gate_nn_x_ = TaylorExpMatrix(this->update_para.tau, ham_nn_x_);
            evolve_gate_nn_y_ = TaylorExpMatrix(this->update_para.tau, ham_nn_y_);
            evolve_gate_nn_z_ = TaylorExpMatrix(this->update_para.tau, ham_nn_z_);
        }

        double SimpleUpdateSweep_(void) override;

        Tensor ham_nn_x_;
        Tensor ham_nn_y_;
        Tensor ham_nn_z_;
        Tensor evolve_gate_H_22_x_;
        Tensor evolve_gate_H_33_x_;
        Tensor evolve_gate_H_23_y_;
        Tensor evolve_gate_H_32_y_;
        Tensor evolve_gate_H_33_y_;
        Tensor evolve_gate_H_22_y_;
        Tensor evolve_gate_H_31_y_;
        Tensor evolve_gate_H_22_z_;
        Tensor evolve_gate_H_33_z_;
        Tensor evolve_gate_H_32_z_;
        Tensor evolve_gate_H_23_z_;
        Tensor evolve_gate_H_31_z_;
        Tensor evolve_gate_H_13_z_;
        Tensor evolve_gate_nn_x_;
        Tensor evolve_gate_nn_y_;
        Tensor evolve_gate_nn_z_;
        double H_;
        double Kz_;
    };


    template<typename TenElemT, typename QNT>
    double KitaevLattice111SimpleUpdateExecutor<TenElemT, QNT>::SimpleUpdateSweep_(void) {
        Timer simple_update_sweep_timer("simple_update_sweep");
        SimpleUpdateTruncatePara para(this->update_para.Dmin, this->update_para.Dmax, this->update_para.Trunc_err);
        double norm = 1.0;
        double e0 = 0.0;

#ifdef QLPEPS_TIMING_MODE
        Timer vertical_nn_projection_x_timer("vertical_nn_projection_x");
#endif
        for (size_t col = 0; col < this->lx_ - 1; col ++){
            if ((col & 1) == 0) {
                for (size_t row = 0; row <= this->ly_ - 1; row = row + 2) {
                    if (row == 0 || row == this->ly_ - 1) {
                        ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_22_x_, {row, col}, HORIZONTAL, para);
                        //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                        norm = result.norm;
                        //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                        e0 += -std::log(norm) / this->update_para.tau;
                    }
                    else {
                        ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_33_x_, {row, col}, HORIZONTAL, para);
                        //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                        norm = result.norm;
                        //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                        e0 += -std::log(norm) / this->update_para.tau;
                    }

                }
            }
            else {
                for (size_t row = 1; row < this->ly_; row = row + 2) {
                    if (row == this->ly_ - 1) {
                        ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_22_x_, {row, col}, HORIZONTAL, para);
                        //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                        norm = result.norm;
                        //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                        e0 += -std::log(norm) / this->update_para.tau;
                    }
                    else {
                        ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_33_x_, {row, col}, HORIZONTAL, para);
                        //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                        norm = result.norm;
                        //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                        e0 += -std::log(norm) / this->update_para.tau;
                    }
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
                    if (col == 0) {
                        if (row == this->ly_-2) {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_22_y_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        } else {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_23_y_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        }
                    } else if (col == this->lx_ - 1) {
                        if (row == this->ly_-2) {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_31_y_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        } else {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_32_y_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        }
                    } else {
                        if (row == this->ly_-2) {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_32_y_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        } else {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_33_y_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        }
                    }

                }
            }
            else {
                for (size_t row = 0; row < this->ly_ - 1; row = row + 2) {
                    if (col == this->lx_ - 1) {
                        if (row == 0) {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_22_y_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        } else if (row == this->ly_-2) {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_31_y_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        } else {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_32_y_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        }
                    } else {
                        if (row == 0) {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_23_y_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        } else if (row == this->ly_-2) {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_32_y_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        } else {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_33_y_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        }
                    }
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
                    if (col == 0) {
                        if (row == 0) {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_22_z_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        } else if (row == this->ly_ - 2) {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_31_z_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        } else {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_32_z_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        }
                    } else if (col == this->lx_-1) {
                        if (row == 0) {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_13_z_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        } else if (row == this->ly_ - 2){
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_22_z_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        } else{
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_23_z_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        }
                    } else {
                        if (row == 0) {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_23_z_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        } else if (row == this->ly_-2) {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_32_z_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        } else {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_33_z_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        }
                    }
                }
            }
            else {
                for (size_t row = 1; row < this->ly_ - 1; row = row + 2) {
                    if (col == this->lx_ -1) {
                        if (row == this->ly_-2) {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_22_z_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        } else {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_23_z_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        }
                    } else {
                        if (row == this->ly_-2) {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_32_z_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        }else {
                            ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_H_33_z_, {row, col}, VERTICAL, para);
                            //ComplexProjectionRes result = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            norm = result.norm;
                            //norm = this->peps_.NearestNeighborSiteProject(evolve_gate_nn_x_, {row, col}, HORIZONTAL, para);
                            e0 += -std::log(norm) / this->update_para.tau;
                        }
                    }
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

#endif //KITAEV_KITAEV_LATTICE_111_SIMPLE_UPDATE_H
