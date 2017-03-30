#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

  VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0) 
  {
		std::cout << "Invalid estimation or ground_truth data" << std::endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i) 
  {
		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse    += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
} // RMSE

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here. Only for Radar measurements.
      Values in Hj are for output var range, direction, range_rate; 
      based on input var px, py, vx, vy
  */
  MatrixXd Hj   = MatrixXd::Zero(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	// temp var 
	float px2 = px * px;
	float py2 = py * py;
  float px2py2 = px2 + py2;
	float pxpy_sqrt = sqrt(px2py2);
	float pxpy_3_2  = pow((px2py2), 1.5);  // AA: Note raising to pow(x, 1.5) gives better result than pow(x, 3/2)

  //check division by zero
  if(fabs(px2 + py2) < 0.00001) {
		std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
		return Hj;
	}

  Hj << px/pxpy_sqrt, py/pxpy_sqrt, 0, 0, 
	      -py/(px2py2), px/(px2py2), 0, 0, 
	      py*(vx*py - vy*px )/pxpy_3_2, px*(vy*px - vx*py)/pxpy_3_2, px/pxpy_sqrt, py/pxpy_sqrt;

 return Hj;
	
} // calcJacobian
