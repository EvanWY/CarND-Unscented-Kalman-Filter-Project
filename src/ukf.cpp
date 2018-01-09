#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  weights = VectorXd(2*n_aug+1);
  weights(0) = lambda / (lambda + n_aug);
  for (int i = 1; i < 2 * n_aug + 1 ; i ++) {
      weights(i) = 0.5 / (lambda + n_aug);
  }

  Xsig_aug_ = MatrixXd(n_aug, 2 * n_aug + 1);
  Xsig_pred_ = MatrixXd(n_x, 2 * n_aug + 1);

  P_.setIdentity();

  cout << "aaaaaaaaa" << P_ << endl;
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    previous_timestamp_ = meas_package.timestamp_;   

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "init with sensor type RADAR" << endl;
      float ro = meas_package.raw_measurements_(0);
      float theta = meas_package.raw_measurements_(1);
      x_ << ro * cos(theta), ro * sin(theta), 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) { 
      cout << "init with sensor type LASER" << endl;
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

	delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = meas_package.timestamp_;

  Prediction(delta_t);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
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
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  AugmentedSigmaPoints();
  SigmaPointPrediction();
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  MatrixXd S = MatrixXd(n_z,n_z);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  PredictRadarMeasurement(Zsig, z_pred, S);
  UpdateRadarState(Zsig, z_pred, S, meas_package.raw_measurements_);
}

void UKF::PredictRadarMeasurement(MatrixXd& Zsig, VectorXd& z_pred, MatrixXd& S) {

  //transform sigma points into measurement space
  for (int i=0; i<2*n_aug+1; i++) {
      double px = Xsig_pred_(0, i);
      double py = Xsig_pred_(1, i);
      double v = Xsig_pred_(2, i);
      double psi = Xsig_pred_(3, i);
      
      Zsig(0, i) = sqrt(px*px+py*py);
      Zsig(1, i) = atan(py / px);
      Zsig(2, i) = (px*cos(psi) + py*sin(psi)) * v / Zsig(0, i);
  }
  //calculate mean predicted measurement
  z_pred = Zsig * weights;
  
  //calculate innovation covariance matrix S
  S << std_radr*std_radr, 0, 0,
        0, std_radphi*std_radphi, 0,
        0, 0, std_radrd*std_radrd;
        
  for (int i=0; i<2*n_aug+1; i++) {
      MatrixXd c = Zsig.col(i) - z_pred;
      S += weights(i) * (c * c.transpose());
  }
}

void UKF::UpdateRadarState(MatrixXd& Zsig, VectorXd& z_pred, MatrixXd& S, VectorXd& z) {
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);
  Tc = MatrixXd::Zero(n_x, n_z);

  //calculate cross correlation matrix
  for (int i=0; i<Zsig.cols(); i++) {
      VectorXd xx = Xsig_pred_.col(i) - x_;
      VectorXd zz = Zsig.col(i) - z_pred;

      while (zz(1)> M_PI) zz(1)-=2.*M_PI;
      while (zz(1)<-M_PI) zz(1)+=2.*M_PI;
      while (xx(3)> M_PI) xx(3)-=2.*M_PI;
      while (xx(3)<-M_PI) xx(3)+=2.*M_PI;

      Tc += weights(i) * xx * zz.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  VectorXd zz = z - z_pred;
  while (zz(1)> M_PI) zz(1)-=2.*M_PI;
  while (zz(1)<-M_PI) zz(1)+=2.*M_PI;

  x_ = x_ + K * zz;
  P_ = P_ - K * S * K.transpose();
}

void UKF::AugmentedSigmaPoints() {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);
 
  //create augmented mean state
  VectorXd aug = VectorXd(n_aug);
  aug << x_ , 0, 0;

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create augmented covariance matrix
  P_aug = P_;
  P_aug.conservativeResize(7, 7);
  P_aug.col(5).setZero();
  P_aug.row(5).setZero();
  P_aug.col(6).setZero();
  P_aug.row(6).setZero();
  P_aug(5, 5) = std_a*std_a;
  P_aug(6, 6) = std_yawdd*std_yawdd;
  
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  
  //create augmented sigma points
  MatrixXd temp_mat = sqrt(lambda + n_aug) * A;
  Xsig_aug_ << aug, temp_mat.colwise() + aug, (-temp_mat).colwise() + aug;
}

void UKF::SigmaPointPrediction() {

  //predict sigma points
  for (int i=0; i<2*n_aug + 1; i++) {
    double dt = delta_t;
    double px = Xsig_aug_(0, i);
    double py = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double psi = Xsig_aug_(3, i);
    double psid = Xsig_aug_(4, i);
    double nu_a = Xsig_aug_(5, i);
    double nu_psidd = Xsig_aug_(6, i);
    
    //write predicted sigma points into right column
    Xsig_pred_(0,i) = px + 0.5*dt*dt*cos(psi)*nu_a
      +( fabs(psid)<0.0001 
          ? v*cos(psi)*dt 
          : v/psid*(sin(psi+psid*dt)-sin(psi)) );
    
    Xsig_pred_(1,i) = py + 0.5*dt*dt*sin(psi)*nu_a
      + ( fabs(psid)<0.0001 
          ? v*sin(psi)*dt 
          : v/psid*(-cos(psi+psid*dt)+cos(psi)) );
    
    Xsig_pred_(2,i) = v+dt*nu_a;
    
    Xsig_pred_(3,i) = psi+psid*dt+0.5*dt*dt*nu_psidd;
    
    Xsig_pred_(4,i) = psid+dt*nu_psidd;
  }
}

void UKF::PredictMeanAndCovariance() {
  //predict state mean
  x_ = Xsig_pred_ * weights;

  //predict state covariance matrix
  MatrixXd XX = Xsig_pred_.colwise() - x_;
  P_ = MatrixXd::Zero(n_x, n_x);
  for (int i=0; i<2*n_aug+1; i++) {
      MatrixXd c = XX.col(i);
      P_ += weights(i) * (c * c.transpose());
  }
}