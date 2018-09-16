/*               
* @Author: Udacity
* @Last Modified by:   debasis123
*/

#include "ukf.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() :
  is_initialized_(false),
  use_laser_(true),
  use_radar_(true),
  n_x_(5),
  n_aug_(n_x_+2),
  lambda_(3.0 - static_cast<double>(n_aug_)),
  n_sigma_(2*n_aug_+1),
  n_z_radar_(3),
  n_z_laser_(2),
  P_(MatrixXd(n_x_, n_x_)),
  Q_(MatrixXd(2, 2)),
  Xsig_pred_(MatrixXd(n_x_, n_sigma_)),
  weights_(VectorXd(n_sigma_)),
  R_radar_(MatrixXd(n_z_radar_, n_z_radar_)),
  R_laser_(MatrixXd(n_z_laser_, n_z_laser_)),
  H_laser_(MatrixXd(2, n_x_)),
  prev_timestamp_(0),
  // assuming the bike stays with an acceleration of [-2, 2] m/s^2, the approx value is 2/1 = 1 m^2/s^4
  // this is based on the suggestion in the lecture
  std_a_(1.0),
  std_yawdd_(0.3),
  // DO NOT MODIFY measurement noise values below, these are provided by the sensor manufacturer.
  std_laspx_(0.15),
  std_laspy_(0.15),
  std_radr_(0.3),
  std_radphi_(0.03),
  std_radrd_(0.3),
  NIS_radar_(0.0),
  NIS_laser_(0.0)
{  
  /**
  TODO:
  Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  */
  
  // create state vector: will be initialized based on the sensor data
  x_ = VectorXd(n_x_);

  // initialize state covariance matrix
  P_.fill(0.0) ;
  // since we assumed the std dev of the lidar x and y measurements to be 0.15
  // this may be initialized more accurately for radar separately
  P_(0,0) = P_(1,1) = 0.15;
  // rest are just assumed correlated
  P_(2,2) = P_(3,3) = P_(4,4) = 1;

  // initialize weights
  const double w0 = lambda_ / (lambda_ + n_aug_);
  const double w1 = 0.5 / (n_aug_ + lambda_);
  weights_(0) = w0;
  for (size_t i = 1; i != n_sigma_; ++i) {
    weights_(i) = w1;
  }

  // initialize predicted sigma point matrix
  Xsig_pred_.fill(0.0);

  // initialize Process Noise Covariance Matrix
  Q_.fill(0.0);
  Q_(0,0) = std_a_ * std_a_;
  Q_(1,1) = std_yawdd_ * std_yawdd_;

  //add measurement noise covariance matrix
  R_radar_.fill(0.0);
  R_radar_(0, 0) = std_radr_ * std_radr_;
  R_radar_(1, 1) = std_radphi_ * std_radphi_;
  R_radar_(2, 2) = std_radrd_ * std_radrd_;

  // add measurement matrix for laser
  H_laser_.fill(0.0);
  H_laser_(0,0) = H_laser_(1,1) = 1;

  // add measurement noise covariance matrix for laser
  R_laser_.fill(0.0);
  R_laser_(0, 0) = std_laspx_ * std_laspx_;
  R_laser_(1, 1) = std_laspy_ * std_laspy_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage& meas_package) {
  /**
  TODO:
  Complete this function! Make sure you switch between lidar and radar measurements.
  */

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    cout << "UKF initialization: " << endl;

    prev_timestamp_ = meas_package.timestamp_;
    
    initialize_UKF(meas_package);

    // done initializing, no need to predict or update
    is_initialized_ = true;

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  
  //compute the time elapsed between the current and previous measurements
  const double dt = (meas_package.timestamp_ - prev_timestamp_) / 1000000.0; //dt - expressed in seconds
  prev_timestamp_ = meas_package.timestamp_;

  // if(dt > 0.0001) {
    Prediction(dt);
  // }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER) {
    UpdateLidar(meas_package);
  }
}


void UKF::initialize_UKF(const MeasurementPackage& meas_package) {
  /**
  TODO:
    * Initialize the state ukf_.x_ with the first measurement.
    * Create the covariance matrix.
    * Remember: you'll need to convert radar from polar to cartesian coordinates.
  */
  // first measurement
  if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER) {

    //set the state with the initial location    
    x_ << meas_package.raw_measurements_(0),
          meas_package.raw_measurements_(1),
          1,
          1,
          0.01;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
    /**
    Convert radar from polar to cartesian coordinates and initialize state.
    */
    const double rho     = meas_package.raw_measurements_(0); // range
    const double phi     = meas_package.raw_measurements_(1); // bearing
    const double rho_dot = meas_package.raw_measurements_(2); // range rate

    /**
    Initialize state.
    */
    x_ << rho * cos(phi),
          rho * sin(phi),
          rho_dot, // not exactly correct value to initialize state due to different frame of ref (Lesson 7, Lec 32)
          phi,
          0.01;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(const double delta_t) {
  /**
  TODO:
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /////////////////////////////////////////////
  // Step 1. Generate Augmented sigma points //
  /////////////////////////////////////////////

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_);

  generate_augmented_sigma_points(Xsig_aug);

  /////////////////////////////////////////
  // Step 2. Predict sigma points at k+1 //
  /////////////////////////////////////////
  
  predict_sigma_points(Xsig_aug, delta_t);

  ///////////////////////////////////////////////
  // Step 3. Compute new mean and sigma points //
  ///////////////////////////////////////////////
  
  compute_mean_covariance_of_sigma_points();
}

void UKF::generate_augmented_sigma_points(MatrixXd& Xsig_aug) {
  // cout << __FUNCTION__ << endl;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug.bottomRightCorner(2, 2) = Q_;

  //create square root matrix
  MatrixXd L = MatrixXd(n_aug_, n_aug_);
  L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (size_t i = 0; i != n_aug_; ++i) { // for each col
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
}

void UKF::predict_sigma_points(const MatrixXd& Xsig_aug, const double delta_t) {
  // cout << __FUNCTION__ << endl;

  // compute for each col states
  for (size_t i = 0; i != n_sigma_; ++i) {
    const double p_x      = Xsig_aug(0, i);
    const double p_y      = Xsig_aug(1, i);
    const double v        = Xsig_aug(2, i);
    const double yaw      = Xsig_aug(3, i);
    const double yawd     = Xsig_aug(4, i);
    const double nu_a     = Xsig_aug(5, i);
    const double nu_yawdd = Xsig_aug(6, i);

    //predicted state values, initialized by values from previous state
    double px_p   = p_x,
           py_p   = p_y,
           v_p    = v,
           yaw_p  = yaw,
           yawd_p = yawd;

    ////////////////////////
    // DETERMINISTIC PART //
    ////////////////////////

    const double ZERO_DIV_LIMIT = 0.001;

    //avoid division by zero
    if (fabs(yawd) > ZERO_DIV_LIMIT) {
      px_p += (v/yawd) * (sin(yaw+yawd*delta_t) - sin(yaw));
      py_p += (v/yawd) * (cos(yaw) - cos(yaw+yawd*delta_t));
    } else {
      px_p += v*delta_t*cos(yaw);
      py_p += v*delta_t*sin(yaw);
    }
    yaw_p += yawd*delta_t;

    /////////////////////
    // STOCHASTIC PART //
    /////////////////////
    
    //add noise
    px_p   += 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p   += 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p    += nu_a * delta_t;
    yaw_p  += 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p += nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }
}

/**
 * angle normalization, to keep the angle between [-pi, pi]
 */
static void normalize_angle(VectorXd& x_diff, const size_t i) {
  while (x_diff(i) > M_PI) {
    x_diff(i) -= 2.0 * M_PI; 
  }
  while (x_diff(i) < -M_PI) {
    x_diff(i) += 2.0 * M_PI; 
  }
}

void UKF::compute_mean_covariance_of_sigma_points() {
  // cout << __FUNCTION__ << endl;

  //predicted state mean
  x_.fill(0.0);
  // go over each sigma point of the predicted sigma points
  for (size_t i = 0; i != n_sigma_; ++i) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  // go over each sigma point of the predicted sigma points
  for (size_t i = 0; i != n_sigma_; ++i) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    normalize_angle(x_diff, 3);

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}


/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage& meas_package) {
  // cout << __FUNCTION__ << endl;

  /**
  TODO:
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the radar NIS.
  */

  ///////////////////////////////////////////////////////////
  // Step 4: Predict Radar measurement mean and covariance //
  ///////////////////////////////////////////////////////////

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_radar_, n_sigma_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar_);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_radar_, n_z_radar_);

  radar_predict_measurement_mean_covariance(Zsig, z_pred, S);

  /////////////////////////////////////////////////
  // Step 5: update current state and covariance //
  /////////////////////////////////////////////////
  
  radar_update_state(meas_package, Zsig, z_pred, S);
}

void UKF::radar_predict_measurement_mean_covariance(Eigen::MatrixXd& Zsig,
                                                    Eigen::VectorXd& z_pred,
                                                    Eigen::MatrixXd& S) {
  // cout << __FUNCTION__ << endl;

  //transform sigma points into measurement space
  for (size_t i = 0; i != n_sigma_; ++i) {
    // extract values for better readibility
    const double p_x = Xsig_pred_(0, i);
    const double p_y = Xsig_pred_(1, i);
    const double v   = Xsig_pred_(2, i);
    const double yaw = Xsig_pred_(3, i);

    const double v1 = cos(yaw) * v;
    const double v2 = sin(yaw) * v;

    // measurement model
    Zsig.fill(0.0);
    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);           //r
    Zsig(1, i) = atan2(p_y, p_x);                   //phi
    Zsig(2, i) = Zsig(0, i) > 0.0001                //r_dot
                 ? (p_x*v1 + p_y*v2 ) / Zsig(0, i)
                 : (p_x*v1 + p_y*v2 ) / 0.0001;
  }

  // predicted measurement mean
  z_pred.fill(0.0);
  for (size_t i = 0; i != n_sigma_; ++i) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // predicted covariance matrix S
  S.fill(0.0);
  for (size_t i = 0; i != n_sigma_; ++i) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization, to keep the angle between [-pi, pi]
    normalize_angle(z_diff, 1);

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  S += R_radar_;
}

void UKF::radar_update_state(const MeasurementPackage& meas_package,
                             const Eigen::MatrixXd& Zsig,
                             const Eigen::VectorXd& z_pred,
                             const Eigen::MatrixXd& S) {
  // cout << __FUNCTION__ << endl;

  //create vector for incoming radar measurement
  VectorXd z = VectorXd(n_z_radar_);
  z << meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);

  //calculate cross-correlation matrix
  Tc.fill(0.0);
  for (size_t i = 0; i != n_sigma_; ++i) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    normalize_angle(z_diff, 1);

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    normalize_angle(x_diff, 3);

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  normalize_angle(z_diff, 1);

  //update state mean and covariance matrix
  x_/*new*/ = x_/*predicted*/ + K * z_diff;
  P_/*new*/ = P_/*predicted*/ - K*S*K.transpose();

  //compute NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage& meas_package) {
  // cout << __FUNCTION__ << endl;

  /**
  TODO:
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the lidar NIS.
  */

  //create vector for incoming lidar measurement
  VectorXd z = VectorXd(2);
  z << meas_package.raw_measurements_;

  // Because lidar measurement update function is a linear function h, 
  // the update step will use the basic Kalman filter equations
  const VectorXd z_pred = H_laser_ * x_;
  const VectorXd z_diff = z - z_pred;
  const MatrixXd Ht     = H_laser_.transpose();
  const MatrixXd S      = H_laser_*P_*Ht + R_laser_;
  const MatrixXd Si     = S.inverse();
  const MatrixXd PHt    = P_ * Ht;
  const MatrixXd K      = PHt * Si;

  //new state
  x_ = x_ + (K*z_diff);

  // new covariance matrix
  const long x_size = x_.size();
  const MatrixXd I  = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K*H_laser_) * P_;

  //compute NIS
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
}

