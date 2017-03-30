#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  // AA - 
  MatrixXd Ft = F_.transpose();
  x_        = F_ * x_;
  P_        = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  // AA:
  VectorXd  z_pred  = H_ * x_;
  VectorXd  y       = z - z_pred;

  MatrixXd  Ht      = H_.transpose();
  MatrixXd  S       = H_ * P_ * Ht + R_;
  MatrixXd  Si      = S.inverse();
  MatrixXd  K       = P_ * Ht * Si;   // Kalman gain
  // new estimates
  x_                = x_ + (K * y);
  long xn           = x_.size();
  MatrixXd  I       = MatrixXd::Identity(xn, xn);
  P_                = (I - K * H_ ) * P_;
} // Update()

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  // AA: calc rho_pred, phi_pred, rhodot_pred
  float rho_pred    = sqrt(pow(x_[0], 2) + pow(x_[1], 2));  
  float phi_pred    = 0.0;
  if (fabs(x_[0]) > 0.001) {
    phi_pred  = atan2(x_[1], x_[0]);    // arctan(py/px)
  }

  float rhodot_pred = 0.0;
  if (fabs(rho_pred) > 0.001) {
    rhodot_pred = (x_[0]*x_[2] + x_[1]*x_[3]) / rho_pred; // (px * vy + py*vx)/rho_pred
  }

  VectorXd  z_pred(3);
  z_pred    << rho_pred, phi_pred, rhodot_pred;

  VectorXd  y   = z - z_pred;     // error in measurement

  MatrixXd  Ht  = H_.transpose();
  MatrixXd  S   = H_ * P_ * Ht + R_;
  MatrixXd  Si  = S.inverse();
  MatrixXd  K   = P_ * Ht * Si;   // Kalman gain

  // new estimates
  x_            = x_ + (K * y);
  long xn       = x_.size();
  MatrixXd  I   = MatrixXd::Identity(xn, xn);
  P_            = (I - K * H_) * P_;

} // UpdateEFK()
