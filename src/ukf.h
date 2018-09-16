/*               
* @Author: Udacity
* @Last Modified by:   debasis123
*/

#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

class UKF {

public:

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(const MeasurementPackage& meas_package);

private:

  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // State dimension
  size_t n_x_;

  // Augmented state dimension
  size_t n_aug_;

  // Sigma point spreading parameter
  double lambda_;

  // number of sigma points
  size_t n_sigma_;

  // measurement dimension, radar can measure r, phi, and r_dot
  size_t n_z_radar_;

  // measurement dimension, laser can measure x, y
  size_t n_z_laser_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // process noise covariance matrix
  Eigen::MatrixXd Q_;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // measurement noise covariance matrix for radar
  Eigen::MatrixXd R_radar_;

  // measurement noise covariance matrix for laser
  Eigen::MatrixXd R_laser_;

  // measurement matrix for laser
  Eigen::MatrixXd H_laser_;

  // time when the state is true, in US
  long long prev_timestamp_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // the current NIS for radar
  double NIS_radar_;

  // the current NIS for laser
  double NIS_laser_;

private:

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(const double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(const MeasurementPackage& meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(const MeasurementPackage& meas_package);


  void initialize_UKF(const MeasurementPackage& meas_package);

  // KF predict
  void generate_augmented_sigma_points(Eigen::MatrixXd& Xsig_aug);
  void predict_sigma_points(const Eigen::MatrixXd& Xsig_aug,
                            const double delta_t);
  void compute_mean_covariance_of_sigma_points();

  // KF update: radar
  void radar_predict_measurement_mean_covariance(Eigen::MatrixXd& Zsig,
                                                 Eigen::VectorXd& z_pred,
                                                 Eigen::MatrixXd& S);
  void radar_update_state(const MeasurementPackage& meas_package,
                          const Eigen::MatrixXd& Zsig,
                          const Eigen::VectorXd& z_pred,
                          const Eigen::MatrixXd& S);
};

#endif /* UKF_H */
