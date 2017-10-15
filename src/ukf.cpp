#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

#define SN 0.001

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  is_initialized_ = false;

  Tools tools;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.25;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma points spreading parameter
  lambda_ = 3 - n_x_;

  // Augmented sigma points spreading parameter
  lambda_aug_ = 3 - n_aug_;

  // Number of sigma points
  n_x_s_ = 2 * n_x_ + 1;

  // Number of augmented sigma points
  n_aug_s_ = 2 * n_aug_ + 1;

  // Number of radar measurement dimension
  n_r_d_ = 3;

  // Pre compute
  n_aug_minus_n_x_ = n_aug_ - n_x_;
  lambda_plus_n_aug_ = lambda_aug_ + n_aug_;
  sqrt_lambda_plus_n_aug_ = sqrt(lambda_plus_n_aug_);
  lambda_division_lambda_plus_n_aug_ = lambda_aug_ / (lambda_plus_n_aug_);
  lambda_plus_n_x_ = lambda_ + n_x_;
  sqrt_lambda_plus_n_x_ = sqrt(lambda_plus_n_x_);

  // Weights of sigma points
  weights_ = VectorXd(n_aug_s_);
  double weight_0 = lambda_division_lambda_plus_n_aug_;
  weights_(0) = weight_0;
  double weight = 0.5 / (lambda_plus_n_aug_);
  for (int i=1; i<n_aug_s_; i++) {
    weights_(i) = weight;
  }

  //TODO: take into consideration measurements covariance
  P_ = MatrixXd::Identity(n_x_, n_x_);

  //measurement matrix
  H_ = MatrixXd(2, n_x_);
  H_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;

  //measurement matrix transpose
  Ht_ = H_.transpose();

  //identidy matrix of 5,5
  I5_ = MatrixXd::Identity(5, 5);

  //laser measurement covariance matrix
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  //radar measurement covariance matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
                      0, std_radphi_ * std_radphi_, 0,
                      0, 0, std_radrd_ * std_radrd_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if ((meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_)
        || (meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_)) {
        return;
  }

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "UKF: " << endl;

    // initialize state vector x
    x_ = VectorXd(n_x_);
    x_ << 0, 0, 0, 0, 0;

    /**
    Initialize state.
    Set the state with the initial position and zero all variables that can't be observed directly.
    */
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian position
      */
      VectorXd ptc(2);
      ptc = tools.ConvertPolarToCartesianPosition(meas_package.raw_measurements_);
      x_ << ptc(0), ptc(1), 0.f, 0.f, 0.f;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0.f, 0.f, 0.f;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  //compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  MatrixXd Xsig = GenerateSigmaPoints();
  MatrixXd Xsig_aug = AugmentedSigmaPoints(Xsig);
  SigmaPointPrediction(Xsig_aug, delta_t);
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  * update the state by using Kalman Filter equations
  */
  VectorXd y = meas_package.raw_measurements_ - H_ * x_;
  MatrixXd S = H_ * P_ * Ht_ + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht_ * Si;

  //new state
  x_ = x_ + (K * y);
  P_ = (I5_ - K * H_) * P_;

  NIS_laser_ = y.transpose() * Si * y;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_r_d_, n_aug_s_);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_r_d_);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_r_d_,n_r_d_);
  //transform sigma points into measurement space
  for (int i=0;i<n_aug_s_;i++) {
      VectorXd x = Xsig_pred_.col(i);

      // extract values for better readability
      double px = x(0);
      double py = x(1);
      double v = x(2);
      double yaw = x(3);

      double rho = sqrt(px * px + py * py);
      double phi = atan2(py,px);
      double rho_dot = (px * cos(yaw) * v + py * sin(yaw) * v) / rho;

      Zsig.col(i) << rho, phi, rho_dot;
  }
  //calculate mean predicted measurement
  z_pred.setZero();
  for (int j=0;j<n_aug_s_;j++) {
      z_pred = z_pred + weights_(j) * Zsig.col(j);
  }
  S.setZero();
  for (int k=0;k<n_aug_s_;k++) {
      VectorXd z_diff = Zsig.col(k) - z_pred;
      z_diff(1) = tools.normalizeAngle(z_diff(1));
      S = S + weights_(k) * z_diff * z_diff.transpose();
  }
  S = S + R_radar_;
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_,n_r_d_);
  Tc.setZero();
  for (int l=0;l<n_aug_s_;l++) {
    //residual
    VectorXd z_diff = Zsig.col(l) - z_pred;
    z_diff(1) = tools.normalizeAngle(z_diff(1));
    // state difference
    VectorXd x_diff = Xsig_pred_.col(l) - x_;
    //angle normalization
    x_diff(3) = tools.normalizeAngle(x_diff(3));

    Tc = Tc + weights_(l) * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  //angle normalization
  z_diff(1) = tools.normalizeAngle(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}

MatrixXd UKF::GenerateSigmaPoints() {
  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, n_x_s_);

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //set first column of sigma point matrix
  Xsig.col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig.col(i+1)      = x_ + sqrt_lambda_plus_n_x_ * A.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt_lambda_plus_n_x_ * A.col(i);
  }

  return Xsig;
}

MatrixXd UKF::AugmentedSigmaPoints(MatrixXd Xsig) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_aug_s_);

  //create augmented mean state
  VectorXd x_noise = VectorXd(2);
  x_noise << 0, 0;
  x_aug.head(n_x_) = x_;
  x_aug.tail(n_aug_ - n_x_) = x_noise;

  //create augmented covariance matrix
  MatrixXd P_noise = MatrixXd(2,2);
  P_noise << std_a_ * std_a_, 0,
            0, std_yawdd_ * std_yawdd_;
  P_aug.setZero();
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(n_aug_minus_n_x_, n_aug_minus_n_x_) = P_noise;

  //create square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();

  //create augmented sigma points
  //set first column of sigma point matrix
  Xsig_aug.col(0) = x_aug;

  //set remaining sigma points
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt_lambda_plus_n_aug_ * A_aug.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt_lambda_plus_n_aug_ * A_aug.col(i);
  }

  return Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd Xsig_aug, double dt) {

  double dt2 = dt * dt;
  double dt2_division_2 = dt2 / 2;

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, n_aug_s_);

  for (int i=0;i < n_aug_s_;i++) {

    // extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // pre compute
    double v_division_yawd = v / yawd;
    double yawd_multiply_dt = yawd * dt;
    double yaw_plus_yawd_multiply_dt = yaw + yawd_multiply_dt;
    double v_multiply_dt = v * dt;

    VectorXd xk = Xsig_aug.col(i);
    VectorXd v1 = xk.head(n_x_);

    VectorXd v2 = VectorXd(n_x_);
    double c1 = fabs(yawd) > SN ? v_division_yawd * (sin(yaw_plus_yawd_multiply_dt) - sin(yaw)) : v_multiply_dt * cos(yaw);
    double c2 = fabs(yawd) > SN ? v_division_yawd * (-cos(yaw_plus_yawd_multiply_dt) + cos(yaw)) : v_multiply_dt * sin(yaw);
    v2 << c1, c2, 0, yawd_multiply_dt, 0;

    VectorXd v3 = VectorXd(n_x_);
    double c3 = dt2_division_2 * cos(yaw) * nu_a;
    double c4 = dt2_division_2 * sin(yaw) * nu_a;
    double c5 = dt * nu_a;
    double c6 = dt2_division_2 * nu_yawdd;
    double c7 = dt * nu_yawdd;
    v3 << c3, c4, c5, c6, c7;

    Xsig_pred.col(i) = v1 + v2 + v3;
  }

  // set sigma point prediction matrix
  Xsig_pred_ = Xsig_pred;
}

void UKF::PredictMeanAndCovariance() {
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);
  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predict state mean
  x.setZero();
  for (int i=0;i<n_aug_s_;i++) {
    x = x + weights_(i) * Xsig_pred_.col(i);
  }
  //predict state covariance matrix
  P.setZero();
  for (int j=0;j<n_aug_s_;j++) {

    //state difference
    VectorXd x_diff = Xsig_pred_.col(j) - x;

    //normalize angle
    x_diff(3) = tools.normalizeAngle(x_diff(3));

    P = P + weights_(j) * x_diff * x_diff.transpose();
  }

  // set predicted state and covariance matrix
  x_ = x;
  P_ = P;
}