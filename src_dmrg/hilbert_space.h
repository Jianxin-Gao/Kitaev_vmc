/*
 * File Name: DefSpinOne.h
 * Description: Define the types, constant and spin-1/2 Hilbert Space.
 * Created by Hao-Xin on 2024/09/21.
 *
 */

#ifndef KITAEV_SRC_DMRG_HILBERT_SPACE_H
#define KITAEV_SRC_DMRG_HILBERT_SPACE_H

#include "qlten/qlten.h"

using TenElemT = qlten::QLTEN_Complex;

using qlten::special_qn::TrivialRepQN;

using qlten::QLTensor;

//namespace kitaev_model {
using QNSctT = qlten::QNSector<TrivialRepQN>;
using IndexT = qlten::Index<TrivialRepQN>;
using Tensor = QLTensor<TenElemT, TrivialRepQN>;
const TrivialRepQN qn0 = TrivialRepQN();
const IndexT pb_out = IndexT({QNSctT(TrivialRepQN(), 2)},
                             qlten::TenIndexDirType::OUT
);
const IndexT pb_in = qlten::InverseIndex(pb_out);

//}


#endif //KITAEV_SRC_DMRG_HILBERT_SPACE_H
