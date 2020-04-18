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

  // newly generated object is uninitilized
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.fill(0.0);
  P_(0,0) = P_(1,1) = P_(2,2) = P_(3,3) = P_(4,4) = 0.5;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.0;    // according to Wikipedia a pro cyclist achieves 2m/s²

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI / 4;  // estimation of personal capabilities - 2PI rotation should take 4s...
  
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

  n_aug_ = 7;
  n_x_ = 5;
  lambda_ = 0;   // was: 3 - n_aug_;

  // create matrices
  weights_ = VectorXd(2*n_aug_+1);
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
}


UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  double dt = (meas_package.timestamp_ - previous_timestamp_)/1000000.0; //dt - expressed in s
  previous_timestamp_ = meas_package.timestamp_;  
  
  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR)
      && (use_radar_==true)){
    /**
       Convert radar from polar to cartesian coordinates and initialize state.
    */
    float px = meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]);
    float py = meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]);
    
    if (!is_initialized_) {
      x_ << px, py, 0, 0, 0; 
      is_initialized_ = true;
      return;
    }

    n_z_ = 3;
    AugmentedSigmaPoints();
    Prediction(dt);
    PredictRadar();
    UpdateRadar(meas_package);
  }
  else if ((meas_package.sensor_type_ == MeasurementPackage::LASER)
	   && (use_laser_==true)){
    /**
       Initialize state.
    */
    if (!is_initialized_) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0; 
      is_initialized_ = true;
      return;
    }

    n_z_ = 2;
    AugmentedSigmaPoints();
    Prediction(dt);
    PredictLidar();
    UpdateLidar(meas_package);    
  }
  return;
}


void UKF::AugmentedSigmaPoints(void) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);
 
  //create sigma point matrix
  // MatrixXd Xsig_aug_ = MatrixXd(n_aug, 2 * n_aug + 1);

  //create augmented mean state
  x_aug = VectorXd::Zero(7);
  x_aug.segment(0, 5) = x_;

  //create augmented covariance matrix
  P_aug = MatrixXd::Zero(7, 7);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug_.col(0) = x_aug;
  for(int i=0; i<n_aug_; ++i)
  {
      Xsig_aug_.col(i+1)        = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
      Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} dt the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
  
  for(int i=0; i<2*n_aug_+1; ++i)
    {
      double p_x      = Xsig_aug_(0, i);
      double p_y      = Xsig_aug_(1, i);
      double v        = Xsig_aug_(2, i);
      double yaw      = Xsig_aug_(3, i);
      double yawd     = Xsig_aug_(4, i);
      double nu_a     = Xsig_aug_(5, i);
      double nu_yawdd = Xsig_aug_(6, i);

      //predict sigma points
      double px_p, py_p;
      if(fabs(yawd)>0.001){
          px_p = p_x + v/yawd * ( sin(yaw + yawd * dt) - sin(yaw));
          py_p = p_y + v/yawd * (-cos(yaw + yawd * dt) + cos(yaw));
      }else{
          //avoid division by zero
          px_p = p_x + v * cos(yaw) * dt;
          py_p = p_y + v * sin(yaw) * dt;
      }
      double v_p = v;
      double yaw_p = yaw + yawd * dt;
      double yawd_p = yawd;
      
      // add noise
      px_p  += 0.5 * nu_a * dt * dt * cos(yaw);
      py_p  += 0.5 * nu_a * dt * dt * sin(yaw);
      v_p   += nu_a * dt;
      yaw_p += 0.5 * nu_yawdd * dt * dt;
      yawd_p += nu_yawdd * dt;

      //write predicted sigma points into right column
      Xsig_pred_(0,i) = px_p;
      Xsig_pred_(1,i) = py_p;
      Xsig_pred_(2,i) = v_p;
      Xsig_pred_(3,i) = yaw_p;
      Xsig_pred_(4,i) = yawd_p;
    }

  //set weights_
  weights_(0) = lambda_ / (lambda_+n_aug_);
  for(int i=1; i<2*n_aug_+1; ++i)
    {
      weights_(i) = 1.0/(2*(lambda_+n_aug_));
    }
  
  //predict state mean
  x_ = Xsig_pred_ * weights_;
  
  //predict state covariance matrix
  P_.fill(0.0); 
  for(int i=0; i<2*n_aug_+1; ++i)
  {
    VectorXd diff = Xsig_pred_.col(i) - x_;
    
    // normalisation of angle
    while(diff(3) >  M_PI) 
      diff(3) -= 2.0*M_PI;
    while(diff(3) < -M_PI)
      diff(3) += 2.0*M_PI;
    
    P_ += weights_(i) * diff * diff.transpose();
  }
}

void UKF::PredictRadar(void){

  //create matrix for sigma points in measurement space
  Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);
  
 //transform sigma points into measurement space
  for(int i=0; i<2*n_aug_ +1; ++i)
  {
    double px     = Xsig_pred_(0, i);
    double py     = Xsig_pred_(1, i);
    double v      = Xsig_pred_(2, i);
    double yaw     = Xsig_pred_(3, i);
    double yaw_dot = Xsig_pred_(4, i);
    
    Zsig_(0, i) = sqrt(px*px + py*py);
    Zsig_(1, i) = atan2(py, px);
    Zsig_(2, i) = (px * cos(yaw) * v + py * sin(yaw) * v) / sqrt(px*px + py*py);
  }

  //calculate mean predicted measurement
  z_pred_ = VectorXd(n_z_);
  z_pred_ = Zsig_ * weights_;

  //calculate innovation covariance matrix S
  S_ = MatrixXd(n_z_, n_z_);
  S_.fill(0.0);
  for(int i=0; i<2*n_aug_+1; ++i)
  {
    VectorXd diff= Zsig_.col(i) - z_pred_;
    //angle normalization
    while (diff(1)> M_PI) diff(1)-=2.*M_PI;
    while (diff(1)<-M_PI) diff(1)+=2.*M_PI;

    S_ += weights_(i) * diff * diff.transpose();
  }

  // add noise matrix R
  R_ = MatrixXd(n_z_, n_z_);
  R_.fill(0.0);
  
  R_(0,0) = std_radr_ * std_radr_;
  R_(1,1) = std_radphi_ * std_radphi_;
  R_(2,2) = std_radrd_ * std_radrd_;
  
  S_ += R_;
}


void UKF::PredictLidar(void){

  // weights are already set in Prediction()
  
  //create matrix for sigma points in measurement space
  Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);
  
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // measurement model
    Zsig_(0,i) = Xsig_pred_(0,i);		//px
    Zsig_(1,i) = Xsig_pred_(1,i);		//py
  }

  //calculate mean predicted measurement
  z_pred_ = VectorXd(n_z_);
  z_pred_ = Zsig_ * weights_;

  //calculate innovation covariance matrix S
  S_ = MatrixXd(n_z_, n_z_);
  S_.fill(0.0);
  for(int i=0; i<2*n_aug_+1; ++i)
  {
    VectorXd diff= Zsig_.col(i) - z_pred_;

    S_ += weights_(i) * diff * diff.transpose();
  }

  // add noise matrix R
  R_ = MatrixXd(n_z_, n_z_);
  R_ << std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;

  S_ += R_;
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  VectorXd z = VectorXd(n_z_);   // Lidar measurement 

  z << meas_package.raw_measurements_[0],   // px
    meas_package.raw_measurements_[1];      // py

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0; i<2*n_aug_+1; ++i)
  {
    VectorXd z_diff = Zsig_.col(i) - z_pred_;    
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S_.inverse();

  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred_;

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S_ * K.transpose();
}


/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  VectorXd z = VectorXd(n_z_);   // radar measurement 

  z << meas_package.raw_measurements_[0], // rho in m
    meas_package.raw_measurements_[1],    // phi in rad
    meas_package.raw_measurements_[2];    // rwho_dot in m/s

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0; i<2*n_aug_+1; ++i)
  {
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

 //calculate Kalman gain K;
  MatrixXd K = Tc * S_.inverse();
  
  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred_;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ +=  K * z_diff;
  P_ -= K * S_ * K.transpose();  
}
