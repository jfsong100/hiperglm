#include <Rcpp.h>
#include <RcppEigen.h>


using Eigen::Map;
using Eigen::VectorXd;
using Eigen::MatrixXd;


// [[Rcpp::export]]
Eigen::VectorXd solve_qr_cpp(Eigen::MatrixXd A, Eigen::VectorXd v) {
  Eigen::VectorXd v_out;
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
  qr.compute(A);
  v_out = qr.solve(v);
  return v_out;
}
