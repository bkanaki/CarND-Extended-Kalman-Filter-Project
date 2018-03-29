#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
  TODO: - Done!
    * predict the state
  */
  // Same for the Laser and Radar as we are not using non-linear prediction
  x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO: - Done!
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO: - Done!
    * update the state by using Extended Kalman Filter equations
  */
  // used for the z_pred
  float h0 = sqrt(x_(0) * x_(0) + x_(1) * x_(1)); 
  float h1 = atan2(x_(1), x_(0));   // docs say angle is limited to [-pi,pi]
  float h2 = (x_(0) * x_(2) + x_(1) * x_(3)) / h0;  // check div(0)?

  std::cout << "Value of angle phi: " << h1 << std::endl;

  VectorXd z_pred = VectorXd(3);
  z_pred << h0, h1, h2;
  VectorXd y = z - z_pred;
  // use the Hj from Jacobian (better to change in process measuremets
  // as need to ensure correct x_ is taken for calculation)
  // MatrixXd H = tools.CalculateJacobian(x_);

  // This is different for laser and radar
  MatrixXd Ht = H_.transpose();

  // rest is the same
  MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

  //new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
