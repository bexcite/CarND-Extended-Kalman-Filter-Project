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
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0,   0,
              0,   0.0009, 0,
              0,   0,   0.09;

  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;


  Hj_ = MatrixXd(3, 4);

  fmt_ = Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " ", "");

  /**
  TODO:
    * Finish initializing the FusionEKF.
  */


  ekf_.x_ = VectorXd(4);
  ekf_.x_ << 0, 0, 0, 0;

  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1000, 0,    0,    0,
             0,    1000, 0,    0,
             0,    0,    1000, 0,
             0,    0,    0,    1000;


  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 0, 0, // 0,2 = dt - will be changed on each step
             0, 1, 0, 0, // 1,3 = dt - will be changed on each step
             0, 0, 1, 0,
             0, 0, 0, 1;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << 0, 0, 0, 0,
             0, 0, 0, 0,
             0, 0, 0, 0,
             0, 0, 0, 0;

  //set the acceleration noise components
  noise_ax = 5;
  noise_ay = 5;


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  cout << ">>>>>>>>>>>>>>>>>\nFEKF::ProcessMeasurement: measurement_pack = " << measurement_pack << endl;

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
    cout << "FEKF::ProcessMeasurement: Not initialized." << endl;


    float x;
    float y;

    if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x = measurement_pack.raw_measurements_[0];
      y = measurement_pack.raw_measurements_[1];

    } else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float ro = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      x = ro * cos(phi);
      y = ro * sin(phi);
    }

    // first measurement
    cout << "FEKF::ProcessMeasurement: EKF Initialization. " << endl;

    // Init state from first measurement
    ekf_.x_ << x, y, 0, 0;

    // Init covariance matrix from first measurement
    ekf_.P_ << 0.3,  0,   0,    0,
               0,    0.3, 0,    0,
               0,    0,   1000, 0,
               0,    0,   0,    1000;



    // timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    /*
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

//      Convert radar from polar to cartesian coordinates and initialize state.

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

//      Initialize state.

    }
  */

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

//  cout << "dt = " << dt << " seconds" <<  endl;



  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
   */

  // Update F according to elapsed time.
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

//  cout << "ekf_.F_ =" << endl << ekf_.F_ << endl;


  // Update the process noise covariance matrix
  //set the process covariance matrix Q
  ekf_.Q_ <<  dt_4/4*noise_ax, 0,               dt_3/2*noise_ax, 0,
              0,               dt_4/4*noise_ay, 0,               dt_3/2*noise_ay,
              dt_3/2*noise_ax, 0,               dt_2*noise_ax,   0,
              0,               dt_3/2*noise_ay, 0,               dt_2*noise_ay;

//  cout << "ekf_.Q_ =" << endl << ekf_.Q_ << endl;



  cout << "FEKF::ProcessMeasurement: Predict " << endl;

  ekf_.Predict();

  cout << "x_ = " << ekf_.x_.format(fmt_) << endl;
  cout << "P_ = " << endl << ekf_.P_ << endl;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  cout << "FEKF::ProcessMeasurement: Update " << endl;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    cout << "Radar Updates === : " << endl;
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    cout << "jacobian ekf_.H_ = " << endl << ekf_.H_ << endl;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    cout << "Laser Updates === :" << endl;
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  // Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " ", "");
  cout << "x_ = " << ekf_.x_.format(fmt_) << endl;
  cout << "P_ = " << endl << ekf_.P_ << endl;
//  cout << "timestamp_ = " << previous_timestamp_ << endl;
  cout << "<<<<<<<<< ProcessMeasurement finished!" << endl;


}
