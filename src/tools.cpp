#include <iostream>
#include "tools.h"
#include <math.h>

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO: - DONE
    * Calculate the RMSE here.
  */

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  size_t N = estimations.size();
  for (size_t i = 0; i < N; ++i) {
    VectorXd d = estimations[i] - ground_truth[i];
    VectorXd dd = d.array() * d.array();
    rmse = rmse + dd;
  }

  rmse = rmse/N;
  rmse = rmse.array().sqrt();

  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO: - DONE
    * Calculate a Jacobian here.
  */

  const double ZERO = 1e-6;
  MatrixXd J(3, 4);

  double px = x_state[0], py = x_state[1];
  double vx = x_state[2], vy = x_state[3];

  double pxpy2 = px*px + py*py;
  if (pxpy2 < ZERO) pxpy2 = ZERO;

  double pxpy2s = sqrt(pxpy2);
  if (pxpy2s < ZERO) pxpy2s = ZERO;

  double pxpy23 = pow(pxpy2, 1.5);
  if (pxpy23 < ZERO) pxpy23 = ZERO;

  double pv = vy*px - vx*py;

  J << px/pxpy2s,     py/pxpy2s,    0,         0,
       -py/pxpy2,     px/pxpy2,     0,         0,
       -py*pv/pxpy23, px*pv/pxpy23, px/pxpy2s, py/pxpy2s;


  return J;
}
