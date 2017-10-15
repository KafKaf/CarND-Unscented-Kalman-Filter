#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to normalize radians to be between -pi and pi
  */
  float normalizeAngle(float angleInRadians);

  /**
  * A helper method to convert polar to cartesian position.
  */
  VectorXd ConvertPolarToCartesianPosition(const VectorXd& x_state);

};

#endif /* TOOLS_H_ */