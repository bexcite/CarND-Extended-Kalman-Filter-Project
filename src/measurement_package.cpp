//
// Created by Pavlo Bashmakov on 3/11/17.
//
#include <iostream>

#include "measurement_package.h"
//#include "Eigen/Dense"

std::ostream& operator<<(std::ostream& os, MeasurementPackage const& mp) {
  std::string sensor_type;
  if (mp.sensor_type_ == MeasurementPackage::LASER) {
    sensor_type = "L";
  } else if (mp.sensor_type_ == MeasurementPackage::RADAR) {
    sensor_type = "R";
  } else {
    sensor_type = "";
  }
  Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " ", "");
  os << sensor_type << mp.raw_measurements_.format(CommaInitFmt) << ", " << mp.timestamp_;
  return os;
}
