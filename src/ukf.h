#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
private:
  Tools tools;
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  // measurement matrix
  Eigen::MatrixXd H_;

  // measurement matrix transpose
  Eigen::MatrixXd Ht_;

  // Identify matrix of 5,5
  Eigen::MatrixXd I5_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Laser measurement noise covariance matrix
  Eigen::MatrixXd R_laser_;

  ///* Radar measurement noise covariance matrix
  Eigen::MatrixXd R_radar_;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma points spreading parameter
  double lambda_;

  ///* Augmented sigma points spreading parameter
  double lambda_aug_;

  ///* Number of sigma points
  int n_x_s_;

  ///* Number of augmented sigma points
  int n_aug_s_;

  ///* Number of radar measurement dimensions
  int n_r_d_;

  ///* NIS calculation for laser
  double NIS_laser_;

  ///* NIS calculation for radar
  double NIS_radar_;

  ///* Pre compute n_aug_ + n_x_
  int n_aug_minus_n_x_;

  ///* Pre compute lambda_aug_ + n_aug_
  int lambda_plus_n_aug_;

  ///* Pre compute sqrt(lambda_aug+ + n_aug_)
  double sqrt_lambda_plus_n_aug_;

  ///* Pre compute lambda / (lambda_ + n_aug_)
  double lambda_division_lambda_plus_n_aug_;

  ///* Pre compute lambda_ + n_x_
  int lambda_plus_n_x_;

  ///* Pre compute sqrt(lambda_ + n_x_)
  double sqrt_lambda_plus_n_x_;

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
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  MatrixXd GenerateSigmaPoints();

  MatrixXd AugmentedSigmaPoints(MatrixXd Xsig);

  void SigmaPointPrediction(MatrixXd Xsig_aug, double dt);

  void PredictMeanAndCovariance();
};

#endif /* UKF_H */
