#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // TODO: YOUR CODE HERE

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  if (0 == estimations.size())
    cout << "Empty estimation vector." << endl;
  //  * the estimation vector size should equal ground truth vector size
  else if (estimations.size() != ground_truth.size())
    cout << "Estimation vector size and ground truth vector size are different." << endl;
  else {
    //accumulate squared residuals
    vector<VectorXd>::const_iterator est = estimations.begin();
    vector<VectorXd>::const_iterator truth = ground_truth.begin();
    do {
      VectorXd    residual = *est - *truth;
      residual = residual.array() * residual.array();
      rmse += residual;
    } while (++est != estimations.end() && ++truth != ground_truth.end());

    //calculate the mean
    rmse /= (double) estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();
  }

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3, 4);
  //recover state parameters
  double  px = x_state(0);
  double  py = x_state(1);
  double  vx = x_state(2);
  double  vy = x_state(3);

  //check division by zero
  double  pxy2 = px * px + py * py;
  if (fabs(pxy2) < 0.0001)
    cout << "Error: both px and py are zero, divide by zero." << endl;
  else {
    //compute the Jacobian matrix
    double  sqrt_pxy2 = sqrt(pxy2);

    double  px_sqrt_pxy2 = px / sqrt_pxy2;
    double  py_sqrt_pxy2 = py / sqrt_pxy2;

    double  pxy2_32 = pxy2 * sqrt_pxy2;

    Hj <<
                            px_sqrt_pxy2,                       py_sqrt_pxy2,            0,            0,
                              -py / pxy2,                          px / pxy2,            0,            0,
      py * (vx * py - vy * px) / pxy2_32, px * (vy * px - vx * py) / pxy2_32, px_sqrt_pxy2, py_sqrt_pxy2;
  }

  return Hj;
}
