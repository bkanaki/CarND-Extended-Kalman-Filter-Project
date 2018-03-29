#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
	rmse << 0,0,0,0;

  // check if the estimations is non empty
  if (estimations.size() < 1) {
	    cout << "Estimations are not present" << endl;
	    return rmse;
	}

  // check if the vector length of estimates and ground truth data match
  if (estimations.size() != ground_truth.size()) {
    cout << "Size of estimations and ground truth do not match" << endl;
    return rmse;
  }

	// accumulate squared residuals
	VectorXd residual(4);
	VectorXd square(4);
	for(int i=0; i < estimations.size(); ++i){
    residual = (estimations[i] - ground_truth[i]);
    square = residual.array() * residual.array();
    rmse = rmse + square;
	}

	// calculate the mean
	rmse = rmse/estimations.size();

	// calculate the squared root
	rmse = rmse.array().sqrt();

	// return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj(3,4);

	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

  // initialize
	// Hj << 0, 0, 0, 0,
	//       0, 0, 0, 0,
	//       0, 0, 0, 0;
	
	float px2py2 = px*px + py*py;
	float sqrt_px2py2 = sqrt(px2py2);
	float vxpy = vx*py;
	float vypx = vy*px;

  //check division by zero error
	if (fabs(px2py2) < 0.0001) {
    cout << "Tools::CalculateJacobian() - Divide by zero error" << endl;
    return Hj;
	}

	//compute the Jacobian matrix
	Hj(0, 0) = px/sqrt_px2py2;
	Hj(0, 1) = py/sqrt_px2py2;
  Hj(0, 2) = 0;
  Hj(0, 3) = 0;
	Hj(1, 0) = -1*py/px2py2;
	Hj(1, 1) = px/px2py2;
  Hj(1, 2) = 0;
  Hj(1, 3) = 0;
	Hj(2, 0) = py*(vxpy - vypx)/(px2py2*sqrt_px2py2);
	Hj(2, 1) = px*(vypx - vxpy)/(px2py2*sqrt_px2py2);
	Hj(2, 2) = px/sqrt_px2py2;
	Hj(2, 3) = py/sqrt_px2py2;

	return Hj;
}
