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
    * Calculate the RMSE here.
  */

	VectorXd rmse(4);
	rmse << 0,0,0,0;
	if (estimations.size() == 0){
		cout << "Tools::CalculateRMSE empty estimations, return." << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
		VectorXd err = estimations[i] - ground_truth[i];
    rmse += VectorXd(err.array() * err.array());
	}

	//calculate the mean
  rmse = rmse * (1.0/estimations.size());

	//calculate the squared root
  rmse = rmse.array().sqrt();
  
  return rmse;
}