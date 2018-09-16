/*               
* @Author: Udacity
* @Last Modified by:   debasis123
*/

#include <uWS/uWS.h>
#include "json.hpp"
#include "ukf.h"
#include "tools.h"

#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;


// for convenience
using json = nlohmann::json;

//////////////////////
// Helper functions //
//////////////////////

/**
 * Check if the SocketIO event has JSON data.
 * If there is data, the JSON object in string format will be returned,
 * else the empty string will be returned.
 */
static std::string hasData(const std::string& s) {
  const auto found_null = s.find("null");
  const auto b1         = s.find_first_of("[");
  const auto b2         = s.find_first_of("]");

  if (found_null != std::string::npos) {
    return "";
  }
  else if (b1 != std::string::npos && b2 != std::string::npos) {
    return s.substr(b1, b2-b1+1);
  }
  return "";
}

/**
 * Gets the LASER data from the simulator
 */
static void fillLaserData(MeasurementPackage& meas_package, istringstream& iss) {
  meas_package.sensor_type_ = MeasurementPackage::SensorType::LASER;
  meas_package.raw_measurements_ = VectorXd(2);
  double px;    iss >> px;
  double py;    iss >> py;
  meas_package.raw_measurements_ << px, py;
  long long timestamp;  iss >> timestamp;
  meas_package.timestamp_ = timestamp;
}

/**
 * Gets RADAR data from the simulator
 */
static void fillRadarData(MeasurementPackage& meas_package, istringstream& iss) {
  meas_package.sensor_type_ = MeasurementPackage::SensorType::RADAR;
  meas_package.raw_measurements_ = VectorXd(3);
  double ro;        iss >> ro;
  double theta;     iss >> theta;
  double ro_dot;    iss >> ro_dot;
  meas_package.raw_measurements_ << ro, theta, ro_dot;
  long long timestamp;  iss >> timestamp;
  meas_package.timestamp_ = timestamp;
}

/**
 * Gets the ground truth values (to be required to compute RMSE)
 */
static void fillGroundTruthValues(VectorXd& gt_values, istringstream& iss) {
  double x_gt;      iss >> x_gt;      gt_values(0) = x_gt;
  double y_gt;      iss >> y_gt;      gt_values(1) = y_gt;
  double vx_gt;     iss >> vx_gt;     gt_values(2) = vx_gt;
  double vy_gt;     iss >> vy_gt;     gt_values(3) = vy_gt;
}

/**
 * Computes estimated values from KF
 */
static void fillEstimatedValues(const UKF& ukf, VectorXd& estimate, istringstream& iss) {
  double p_x = ukf.x_(0);       estimate(0) = p_x;
  double p_y = ukf.x_(1);       estimate(1) = p_y;

  double v   = ukf.x_(2);
  double yaw = ukf.x_(3);

  double v1 = cos(yaw) * v;
  double v2 = sin(yaw) * v;
  
  estimate(2) = v1;
  estimate(3) = v2;
}

/**
 * Computes data to be fedback to the simulator from C++
 */
static void fillFeedbackToSimulator(json& msgJson, const VectorXd& estimate, const VectorXd& RMSE) {
  // kalman filter estimated position x and y
  // msgJson["estimate_x"] = p_x;
  msgJson["estimate_x"] = estimate(0);
  // msgJson["estimate_y"] = p_y;
  msgJson["estimate_y"] = estimate(1);
  msgJson["rmse_x"]     = RMSE(0);
  msgJson["rmse_y"]     = RMSE(1);
  msgJson["rmse_vx"]    = RMSE(2);
  msgJson["rmse_vy"]    = RMSE(3);
}


/////////////////////////////////////
// Main entry point of the project //
/////////////////////////////////////

/**
 * communicates with the Term 2 Simulator receiving data measurements,
 * calls a function to run the Kalman filter,
 * calls a function to calculate RMSE
 */
int main()
{
  uWS::Hub h;

  // Create a Kalman Filter instance
  UKF ukf;

  // used to compute the RMSE later
  Tools tools;
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;

  ///////////////////////////////////////////////////////
  // BELOW WE DEFINE A SERIES OF FUNCTION OBJECTS on h //
  ///////////////////////////////////////////////////////

  h.onMessage(
    [&ukf, &tools, &estimations, &ground_truth] (uWS::WebSocket<uWS::SERVER> ws,
                                                 char *data,
                                                 size_t length,
                                                 uWS::OpCode opCode) {
    ////////////////////////////////////////////////////////
    // Check if there is incoming data from the simulator //
    // Read the simulation data and parse, if any         //
    ////////////////////////////////////////////////////////

    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event

    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(std::string(data));

      if (!s.empty()) {
        auto j = json::parse(s);
        std::string event = j[0].get<std::string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          // INPUT: values provided by the simulator to the c++ program
          // the measurement that the simulator observed (either lidar or radar)
          string sensor_measurment = j[1]["sensor_measurement"];
          
          MeasurementPackage meas_package;
          istringstream iss(sensor_measurment);

          // read first element from the current line to determine whether it's a Lidar or Radar data 
      	  string sensor_type;
      	  iss >> sensor_type;

      	  if (sensor_type.compare("L") == 0) {
            fillLaserData(meas_package, iss);
          }
          else if (sensor_type.compare("R") == 0) {
            fillRadarData(meas_package, iss);
          }

          //Call ProcessMeasurment(meas_package) for Kalman filter
          ukf.ProcessMeasurement(meas_package);       


          //////////////////////////////////////////////////
          // Compute the RMSE values and Estimated values //
          //////////////////////////////////////////////////


          // ground truth values
          VectorXd gt_values(4);
          fillGroundTruthValues(gt_values, iss);            
          ground_truth.push_back(gt_values);
          
          // Push the current estimated x,y positon from the Kalman filter's state vector
          VectorXd estimate(4);
          fillEstimatedValues(ukf, estimate, iss);
          estimations.push_back(estimate);

          // compute RMSE for comparison
      	  VectorXd RMSE = tools.CalculateRMSE(estimations, ground_truth);


          /////////////////////////////////////////////////////////////////
          // OUTPUT: values provided by the c++ program to the simulator //
          /////////////////////////////////////////////////////////////////
          json msgJson;
          fillFeedbackToSimulator(msgJson, estimate, RMSE);
          std::string msg = "42[\"estimate_marker\"," + msgJson.dump() + "]";
          // std::cout << msg << std::endl;
          // std:: cout << "RMSE values: ";
          // std::cout << "[p_x: " << RMSE(0) << " p_y: " << RMSE(1) 
          //           << " vx: " << RMSE(2) << " vy: " << RMSE(3) << "]" 
          //           << std::endl;
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      }
      // o/w, the SocketIO event has no JSON data.
      else {
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  }); // end function object

  // We don't need this since we're not using HTTP but if it's removed the program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }

  h.run();
}
