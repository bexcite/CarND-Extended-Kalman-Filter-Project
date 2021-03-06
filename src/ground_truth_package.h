#ifndef GROUND_TRUTH_PACKAGE_H_
#define GROUND_TRUTH_PACKAGE_H_

#include "Eigen/Dense"

class GroundTruthPackage {
public:
  long timestamp_;

  enum SensorType{
    LASER,
    RADAR
  } sensor_type_;

  Eigen::VectorXd gt_values_;

};


std::ostream& operator<<(std::ostream& os, GroundTruthPackage const& gt);

#endif /* GROUND_TRUTH_PACKAGE_H_ */
