#include <iostream>
#include "kalman_filter.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;  // state(location, velocity)
  P_ = P_in;  // object covariance matrix (uncertainty)
  F_ = F_in;  // state transition matrix
  H_ = H_in;  // measurement matrix
  R_ = R_in;  // measurement covariance matrix (uncertainty)
  Q_ = Q_in;  // process covariance matrix (uncertainty)
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  // Calculate error term, y
  VectorXd  z_pred = H_ * x_;
  VectorXd  y = z - z_pred;

  UpdateHelper(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  double  p_x = x_[0];
  double  p_y = x_[1];
  double  v_x = x_[2];
  double  v_y = x_[3];

  double  rho = sqrt(p_x * p_x + p_y * p_y);
  double  phi = atan2(p_y, p_x);
  if (rho < 0.0001) // Avoid division by zero.
    rho = 0.0001;
  double  rho_dot = (p_x * v_x + p_y * v_y) / rho;

  VectorXd z_pred = VectorXd(3);  // h(x)
  z_pred << rho, phi, rho_dot;

  // Error term, y
  VectorXd y = z - z_pred;

  // Normalize phi to the range [-pi, pi]
  if (y[1] > 0.0)
    while (y[1] > M_PI) y[1] -= (2 * M_PI);
  else
    while (y[1] < -M_PI) y[1] += (2 * M_PI);

  UpdateHelper(y);
}

void KalmanFilter::UpdateHelper(const VectorXd &y) {

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;  // Kalman gain

  // new estimate
  x_ = x_ + (K * y);
  long  x_size = (long) x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - (K * H_)) * P_;
}