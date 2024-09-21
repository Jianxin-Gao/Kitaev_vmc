
#include "hilbert_space.h"

Tensor sigma_x = Tensor({pb_in, pb_out});
Tensor sigma_y = Tensor({pb_in, pb_out});
Tensor sigma_z = Tensor({pb_in, pb_out});
Tensor id = Tensor({pb_in, pb_out});

void OperatorInitial() {
  static bool initialized = false;
  if (!initialized) {
    sigma_z({0, 0}) = 1.0;
    sigma_z({1, 1}) = -1.0;

    sigma_x({1, 0}) = 1.0;
    sigma_x({0, 1}) = 1.0;

    sigma_y({0, 1}) = std::complex<double>(0, 1);
    sigma_y({1, 0}) = -std::complex<double>(0, 1);
    //(0, -i; i ,0)

    id({0, 0}) = 1;
    id({1, 1}) = 1;

    initialized = true;
  }
}