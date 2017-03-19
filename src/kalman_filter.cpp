#include "kalman_filter.h"
#include <math.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO: - DONE
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO: - DONE
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

VectorXd KalmanFilter::CartToPolar(const VectorXd &x) {
  VectorXd pol(3);

  pol[0] = sqrt(x[0] * x[0] + x[1] * x[1]);

  pol[1] = atan2(x[1], x[0]);

//  std::cout << "atan2 = " << pol[1] << std::endl;

  double d = pol[0];
  if (d < 1e-6) d = 1e-6;

  pol[2] = (x[0] * x[2] + x[1] * x[3]) / d;
  return pol;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO: - DONE
    * update the state by using Extended Kalman Filter equations
  */

  //VectorXd z_pred = H_ * x_;
  VectorXd z_pred = CartToPolar(x_);


  VectorXd y = z - z_pred;

  // Make it in -pi;pi range
  // but there actually no such data in provided datasets :(
  while (std::abs(y[1]) > M_PI) {
    if (y[1] < 0) y[1] += M_PI;
    else if (y[1] > 0) y[1] -= M_PI;
    //std::cout << "y_phi = " << y[1] << " pi" << M_PI << std::endl;
  }

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}
