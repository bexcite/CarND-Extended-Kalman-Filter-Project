//
// Created by Pavlo Bashmakov on 3/11/17.
//

#include "ground_truth_package.h"

std::ostream& operator<<(std::ostream& os, GroundTruthPackage const& gt) {
//  std::string sensor_type;
//  if (gt.sensor_type_ == GroundTruthPackage::LASER) {
//    sensor_type = "L";
//  } else if (gt.sensor_type_ == GroundTruthPackage::RADAR) {
//    sensor_type = "R";
//  } else {
//    sensor_type = "";
//  }
  Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " ", "");
  os /*<< sensor_type*/ << gt.gt_values_.format(CommaInitFmt);
  return os;
}
