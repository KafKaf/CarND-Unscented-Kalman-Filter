#include <iostream>
#include "tools.h"

#define PI 3.14159265

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE
  */

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size()
          || estimations.size() == 0){
      cout << "Invalid estimation or ground_truth data" << endl;
      return rmse;
  }

  //accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){

      VectorXd residual = estimations[i] - ground_truth[i];

      //coefficient-wise multiplication
      residual = residual.array()*residual.array();
      rmse += residual;
  }

  //calculate the mean
  rmse = rmse/estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

float Tools::normalizeAngle(float angleInRadians) {
  while (angleInRadians < -PI) {
    angleInRadians += 2.f*PI;
  };

  while (angleInRadians > PI) {
    angleInRadians -= 2.f*PI;
  }

  return angleInRadians;
}

VectorXd Tools::ConvertPolarToCartesianPosition(const VectorXd& x_state) {
  //recover state parameters
  float rho = x_state(0);
  float phi = x_state(1);

  // find x,y position from ro and phi
  float px = rho * cos(phi);
  float py = rho * sin(phi);

  VectorXd ptc(2);
  ptc << px, py;

  return ptc;
}