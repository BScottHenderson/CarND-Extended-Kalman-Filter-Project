#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  //measurement covariance matrix - laser
  R_laser_ <<
    0.0225, 0,
    0,      0.0225;

  //measurement covariance matrix - radar
  R_radar_ <<
    0.09, 0,      0,
    0,    0.0009, 0,
    0,    0,      0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  H_laser_ <<
    1, 0, 0, 0,
    0, 1, 0, 0;

  /*
  Initialize the KalmanFilter object.
  */

  // Initial state.
  VectorXd  x = VectorXd(4);
  x << 0, 0, 0, 0;

  // State covariance matrix P
  MatrixXd  P = MatrixXd(4, 4);
  P <<
    1, 0,    0,    0,
    0, 1,    0,    0,
    0, 0, 1000,    0,
    0, 0,    0, 1000;

  // Initial transition matrix F
  MatrixXd  F = MatrixXd(4, 4);
  float     t = 1.0;  // The time value, t, will be updated for each measurement.
  F <<
    1, 0, t, 0,
    0, 1, 0, t,
    0, 0, 1, 0,
    0, 0, 0, 1;

  // The measurement matrix and covariance matrix are reset for each measurement
  // based on measurement type so don't bother setting them here.
  MatrixXd  H;  // Measurement matrix.
  MatrixXd  R;  // Measurement covariance matrix (uncertainty).

  // Process covariance matrix, Q. This will be updated for each measurement so
  // just set it to zero for now.
  MatrixXd  Q = MatrixXd(4, 4);
  Q <<
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0;

  ekf_.Init(x, P, F, H, R, Q);

  // Acceleration noise:
  noise_ax_ = 9.0;  // /sigma^2_{ax}
  noise_ay_ = 9.0;  // /sigma^2_{ay}
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */

    // first measurement
    cout << "EKF: " << endl;

    // Save the current timestamp so we can calculate delta next time.
    previous_timestamp_ = measurement_pack.timestamp_;
    current_step_ = 0;

    // Initial state (location, velocity).
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      // Measurement vector: (rho, phi, rho_dot)
      double  rho = measurement_pack.raw_measurements_[0];
      double  phi = measurement_pack.raw_measurements_[1];
      double  rho_dot = measurement_pack.raw_measurements_[2];

      // State vector: (p_x, p_y, v_x, v_y) -> same as for laser measurement
      if (phi > 0.0)  // Normalize phi to the range [-pi, pi]
        while (phi > M_PI) phi -= (2 * M_PI);
      else
        while (phi < -M_PI) phi += (2 * M_PI);
      double  cos_phi = cos(phi);
      double  sin_phi = sin(phi);

      double  p_x = rho * cos_phi;      // Translate polar coordinates to Cartesian.
      double  p_y = rho * sin_phi;
      double  v_x = 0; // rho_dot * cos_phi;
      double  v_y = 0; // rho_dot * sin_phi;

      ekf_.x_ << p_x, p_y, v_x, v_y;
      cout << "Initial measurement: RADAR" << endl;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */

      // Set the state with the initial location and zero velocity.
      double  p_x = measurement_pack.raw_measurements_[0];
      double  p_y = measurement_pack.raw_measurements_[1];
      double  v_x = 0;
      double  v_y = 0;

      ekf_.x_ << p_x, p_y, v_x, v_y;
      cout << "Initial measurement: LASER" << endl;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  //cout << "step: " << ++current_step_ << endl;

  // Calculate elapsed time in seconds.
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0F;
  previous_timestamp_ = measurement_pack.timestamp_;

  // Modify the state transition matrix, F, so that the time is integrated.
  ekf_.F_(0, 2) = ekf_.F_(1, 3) = dt;

  // Set the process covariance matrix, Q, based on elapsed time and noise.
  float   dt2 = dt * dt;
  float   dt3 = dt2 * dt;
  float   dt4 = dt3 * dt;
  ekf_.Q_ <<
    dt4 / 4 * noise_ax_,                   0, dt3 / 2 * noise_ax_,                   0,
                      0, dt4 / 4 * noise_ay_,                   0, dt3 / 2 * noise_ay_,
    dt3 / 2 * noise_ax_,                   0,     dt2 * noise_ax_,                   0,
                      0, dt3 / 2 * noise_ay_,                   0,     dt2 * noise_ay_;

  // Call the Kalman Filter Predict() function.
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  //cout << "x_ = " << endl << ekf_.x_ << endl;
  //cout << "P_ = " << endl << ekf_.P_ << endl;
}
