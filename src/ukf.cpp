#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 5.0 / 10.0; // TODO: fixme

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI / 20.0; // TODO: fixme

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;

  n_x_ = 5;

  n_aug_ = n_x_ + 2;

  lambda_ = 3 - n_aug_;

  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  P_ = MatrixXd(n_x_, n_x_);
  P_.fill(0.0);

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  time_us_ = 0;

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.0);

  NIS_radar_ = 0.0;
  NIS_laser_ = 0.0;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_) {

    x_.fill(0.0);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      x_ << meas_package.raw_measurements_[0] * std::cos(meas_package.raw_measurements_[1]),
          meas_package.raw_measurements_[0] * std::sin(meas_package.raw_measurements_[1]),
          0,
          0,
          0;

    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << meas_package.raw_measurements_[0],
          meas_package.raw_measurements_[1],
          0,
          0,
          0;
    }

    // Initialize covariance matrix
    P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, M_PI_2, 0,
        0, 0, 0, 0, 5;

    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;
    return;
  }

  const double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;    //delta_t - expressed in seconds
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  } else {
    std::cout << "Measurement type is not defined!" << std::endl;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);

  x_aug.head(x_.size()) = x_;

  P_aug.topLeftCorner(P_.rows(), P_.cols()) = P_;

  MatrixXd Q = Eigen::MatrixXd::Zero(2, 2);
  Q(0, 0) = std_a_ * std_a_;
  Q(1, 1) = std_yawdd_ * std_yawdd_;

  P_aug.bottomRightCorner(Q.rows(), Q.cols()) = Q;


  //std::cout << "P_aug: " << std::endl << P_aug << std::endl;

  MatrixXd A = P_aug.llt().matrixL();

  const double sqrtLambdaNx = std::sqrt(lambda_ + n_aug_);
  int j = 0;
  Xsig_aug.col(j++) = x_aug;

  for (int k = 0; j <= n_aug_; ++j, ++k) {
    Xsig_aug.col(j) = x_aug + A.col(k) * sqrtLambdaNx;
    Xsig_aug.col(j) = x_aug + A.col(k) * sqrtLambdaNx;
  }

  for (int k = 0; j <= 2 * n_aug_; ++j, ++k) {
    Xsig_aug.col(j) = x_aug - A.col(k) * sqrtLambdaNx;
    Xsig_aug.col(j) = x_aug - A.col(k) * sqrtLambdaNx;
  }

  //predict sigma points
  auto PredictSigmaPoints = [this, &delta_t](const VectorXd &xk, const bool &isNonZero) -> VectorXd {
    VectorXd xd = VectorXd(n_x_);
    xd.fill(0.0);

    const double &x = xk(0);
    const double &y = xk(1);
    const double &v = xk(2);
    const double &phi = xk(3);
    const double &phi_dot = xk(4);
    const double &v_a = xk(5);
    const double &v_phi_dot_dot = xk(6);

    const double sin_phi = sin(phi);
    const double cos_phi = cos(phi);
    const double phi_dot_delta_t = phi_dot * delta_t;

    const double half_delta_t_square = delta_t * delta_t / 2.0;

    xd(0) = (isNonZero ? (v / phi_dot) * (+sin(phi + phi_dot_delta_t) - sin_phi) : v * cos_phi * delta_t)
        + (half_delta_t_square * cos_phi * v_a);
    xd(1) = (isNonZero ? (v / phi_dot) * (-cos(phi + phi_dot_delta_t) + cos_phi) : v * sin_phi * delta_t)
        + (half_delta_t_square * sin_phi * v_a);
    xd(2) = delta_t * v_a;
    xd(3) = phi_dot_delta_t + half_delta_t_square * v_phi_dot_dot;
    xd(4) = delta_t * v_phi_dot_dot;

    return xk.head(xd.rows()) + xd;
  };

  Xsig_pred_.fill(0.0);

  for (int j = 0; j < Xsig_pred_.cols(); ++j) {
    const VectorXd &xk = Xsig_aug.col(j);
    const bool isNonZero = fabs(xk(4)) > 1e-4;
    VectorXd xkp1 = PredictSigmaPoints(xk, isNonZero);
    xkp1(3) = AngleNorm(xkp1(3));
    Xsig_pred_.col(j) = xkp1;
  }

  //set weights
  weights_.fill(1.0);
  weights_ *= 0.5 / (lambda_ + n_aug_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);


  //predict state mean
  x_ = Xsig_pred_ * weights_;

  x_(3) = AngleNorm(x_(3));

//predict state covariance matrix
  MatrixXd Xsig_pred_minus_x = Xsig_pred_.colwise() - x_;

  for (int j = 0; j < Xsig_pred_minus_x.cols(); ++j) {
    Xsig_pred_minus_x(3, j) = AngleNorm(Xsig_pred_minus_x(3, j));
  }

  P_ = Xsig_pred_minus_x * weights_.asDiagonal() * Xsig_pred_minus_x.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  //set measurement dimension, Laser can measure pos1, and pos2
  const int n_z = 2;

  auto TransformToMeasurementSpace = [&n_z](const VectorXd &xk) {
    VectorXd z_pred = VectorXd(n_z);
    z_pred(0) = xk(0);
    z_pred(1) = xk(1);
    return z_pred;
  };

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  for (int j = 0; j < Xsig_pred_.cols(); ++j) {
    Zsig.col(j) = TransformToMeasurementSpace(Xsig_pred_.col(j));
  }

  //predict state mean
  z_pred = Zsig * weights_;

  const MatrixXd Zsig_minus_z_pred = Zsig.colwise() - z_pred;

  MatrixXd R(z_pred.rows(), z_pred.rows());
  R.fill(0.0);
  R(0, 0) = std_laspx_ * std_laspx_;
  R(1, 1) = std_laspy_ * std_laspy_;

  S = Zsig_minus_z_pred * weights_.asDiagonal() * Zsig_minus_z_pred.transpose() + R;

  NIS_laser_ = UpdateState(Zsig, z_pred, S, meas_package.raw_measurements_);

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //set measurement dimension, radar can measure r, phi, and r_dot
  const int n_z = 3;

  auto TransformToMeasurementSpace = [&n_z](const VectorXd &xk) {
    const double &px = xk(0);
    const double &py = xk(1);
    const double &v = xk(2);
    const double &phi = xk(3);

    double c1 = px * px + py * py;
    if (fabs(c1) < 1e-4) {
      c1 = 1e-4; // some small number and positive
    }
    const double c2 = sqrt(c1);

    VectorXd z_pred = VectorXd(n_z);
    z_pred(0) = c2;
    z_pred(1) = std::atan2(py, px);

    // stability check
    z_pred(2) = (px * cos(phi) * v + py * sin(phi) * v) / c2;

    return z_pred;
  };

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  for (int j = 0; j < Xsig_pred_.cols(); ++j) {
    Zsig.col(j) = TransformToMeasurementSpace(Xsig_pred_.col(j));
  }

  //predict state mean
  z_pred = Zsig * weights_;

  MatrixXd Zsig_minus_z_pred = Zsig.colwise() - z_pred;

  for (int j = 0; j < Zsig_minus_z_pred.cols(); ++j) {
    Zsig_minus_z_pred(1, j) = AngleNorm(Zsig_minus_z_pred(1, j));
  }

  MatrixXd R(z_pred.rows(), z_pred.rows());
  R.fill(0.0);
  R(0, 0) = std_radr_ * std_radr_;
  R(1, 1) = std_radphi_ * std_radphi_;
  R(2, 2) = std_radrd_ * std_radrd_;

  S = Zsig_minus_z_pred * weights_.asDiagonal() * Zsig_minus_z_pred.transpose() + R;

  NIS_radar_ = UpdateState(Zsig, z_pred, S, meas_package.raw_measurements_);

}

// angle normalization
double UKF::AngleNorm(double x) {
  while (x > +M_PI) x -= 2. * M_PI;
  while (x < -M_PI) x += 2. * M_PI;
  return x;
};

double UKF::UpdateState(const MatrixXd &Zsig,
                        const VectorXd &z_pred,
                        const MatrixXd &S,
                        const VectorXd &z) {

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, z.rows());

  MatrixXd Xsig_pred_minus_x = Xsig_pred_.colwise() - x_;

  for (int j = 0; j < Xsig_pred_minus_x.cols(); ++j) {
    Xsig_pred_minus_x(3, j) = AngleNorm(Xsig_pred_minus_x(3, j));
  }

  MatrixXd Zsig_minus_z_pred = Zsig.colwise() - z_pred;

  if (Zsig.rows() == 3) {
    for (int j = 0; j < Zsig_minus_z_pred.cols(); ++j) {
      Zsig_minus_z_pred(1, j) = AngleNorm(Zsig_minus_z_pred(1, j));
    }
  }

  Tc = Xsig_pred_minus_x * weights_.asDiagonal() * Zsig_minus_z_pred.transpose();

  const MatrixXd S_inv = S.inverse();
  const MatrixXd K = Tc * S_inv;

  x_ += K * (z - z_pred);
  P_ -= K * S * K.transpose();

  // Normalized Innovation Squared (NIS)
  const VectorXd z_pred_minus_z = z_pred - z;
  return z_pred_minus_z.transpose() * S_inv * z_pred_minus_z;
}
