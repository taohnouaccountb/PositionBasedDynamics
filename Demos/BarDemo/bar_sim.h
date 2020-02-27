#include "Common/Common.h"
// #include <vector>
// #include <Eigen/Dense>
// typedef double Real;
// using Vector3r = Eigen::Matrix<Real, 3, 1, Eigen::DontAlign>;

void sim_init(int argc, char **argv);
Vector3r sim_exec(const std::vector<Eigen::Isometry3d> *transforms, const Eigen::VectorXd *dt_inv);
void sim_test();