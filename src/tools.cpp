/*               
* @Author: Udacity
* @Last Modified by:   debasis123
*/

#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd>& estimations,
                              const vector<VectorXd>& ground_truth) {
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  eigen_assert(estimations.size() != 0 
    && "estimations vector size cannot be 0");
  eigen_assert(estimations.size() == ground_truth.size() 
    && "estimations and ground_truth size should be same");

  //accumulate squared residuals
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  for (size_t i=0; i != estimations.size(); ++i) {
    VectorXd residual = estimations.at(i) - ground_truth.at(i);
    //coefficient-wise multiplication
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse /= estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}
