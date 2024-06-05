

#ifndef SPIN_COMPLEX_H
#define SPIN_COMPLEX_H

#include "qlten/qlten.h"

using TenElemT = qlten::QLTEN_Complex;
using qlten::special_qn::U1QN;
using qlten::QLTensor;

using QNSctT = qlten::QNSector<U1QN>;
using IndexT = qlten::Index<U1QN>;
using Tensor = QLTensor<TenElemT, U1QN>;
const U1QN qn0 = U1QN(0); //N(particle number), Sz

const IndexT pb_out = IndexT({
                                 QNSctT(U1QN(0), 2)},
                             qlten::TenIndexDirType::OUT
);

const IndexT pb_in = qlten::InverseIndex(pb_out);

std::string peps_path = "peps";

#endif //SPIN_COMPLEX_H
