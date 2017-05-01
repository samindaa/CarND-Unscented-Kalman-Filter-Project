
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "Eigen/Dense"
#include "ukf.h"
#include "ground_truth_package.h"
#include "measurement_package.h"
#include "cmaes/cmaes_interface.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

void check_arguments(int argc, char *argv[]) {
  string usage_instructions = "Usage instructions: ";
  usage_instructions += argv[0];
  usage_instructions += " path/to/input.txt output.txt";

  bool has_valid_args = false;

  // make sure the user has provided input and output files
  if (argc == 1) {
    cerr << usage_instructions << endl;
  } else if (argc == 2) {
    cerr << "Please include an output file.\n" << usage_instructions << endl;
  } else if (argc == 3) {
    has_valid_args = true;
  } else if (argc > 3) {
    cerr << "Too many arguments.\n" << usage_instructions << endl;
  }

  if (!has_valid_args) {
    exit(EXIT_FAILURE);
  }
}

void check_files(ifstream &in_file, string &in_name,
                 ofstream &out_file, string &out_name) {
  if (!in_file.is_open()) {
    cerr << "Cannot open input file: " << in_name << endl;
    exit(EXIT_FAILURE);
  }

  if (!out_file.is_open()) {
    cerr << "Cannot open output file: " << out_name << endl;
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char *argv[]) {

  check_arguments(argc, argv);

  string in_file_name_ = argv[1];
  ifstream in_file_(in_file_name_.c_str(), ifstream::in);

  string out_file_name_ = argv[2];
  ofstream out_file_(out_file_name_.c_str(), ofstream::out);

  check_files(in_file_, in_file_name_, out_file_, out_file_name_);

  /**********************************************
   *  Set Measurements                          *
   **********************************************/

  vector<MeasurementPackage> measurement_pack_list;
  vector<GroundTruthPackage> gt_pack_list;

  string line;

  // prep the measurement packages (each line represents a measurement at a
  // timestamp)
  while (getline(in_file_, line)) {
    string sensor_type;
    MeasurementPackage meas_package;
    GroundTruthPackage gt_package;
    istringstream iss(line);
    long long timestamp;

    // reads first element from the current line
    iss >> sensor_type;

    if (sensor_type.compare("L") == 0) {
      // laser measurement

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::LASER;
      meas_package.raw_measurements_ = VectorXd(2);
      float px;
      float py;
      iss >> px;
      iss >> py;
      meas_package.raw_measurements_ << px, py;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    } else if (sensor_type.compare("R") == 0) {
      // radar measurement

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);
      float ro;
      float phi;
      float ro_dot;
      iss >> ro;
      iss >> phi;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro, phi, ro_dot;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    }

    // read ground truth data to compare later
    float x_gt;
    float y_gt;
    float vx_gt;
    float vy_gt;
    iss >> x_gt;
    iss >> y_gt;
    iss >> vx_gt;
    iss >> vy_gt;
    gt_package.gt_values_ = VectorXd(4);
    gt_package.gt_values_ << x_gt, y_gt, vx_gt, vy_gt;
    gt_pack_list.push_back(gt_package);
  }

  // Create a UKF instance
  UKF ukf;

  // used to compute the RMSE later
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;

  // start filtering from the second frame (the speed is unknown in the first
  // frame)

  size_t number_of_measurements = measurement_pack_list.size();

  // column names for output file
  out_file_ << "time_stamp" << "\t";
  out_file_ << "px_state" << "\t";
  out_file_ << "py_state" << "\t";
  out_file_ << "v_state" << "\t";
  out_file_ << "yaw_angle_state" << "\t";
  out_file_ << "yaw_rate_state" << "\t";
  out_file_ << "sensor_type" << "\t";
  out_file_ << "NIS" << "\t";
  out_file_ << "px_measured" << "\t";
  out_file_ << "py_measured" << "\t";
  out_file_ << "px_ground_truth" << "\t";
  out_file_ << "py_ground_truth" << "\t";
  out_file_ << "vx_ground_truth" << "\t";
  out_file_ << "vy_ground_truth" << "\n";

  for (size_t k = 0; k < number_of_measurements; ++k) {
    // Call the UKF-based fusion
    ukf.ProcessMeasurement(measurement_pack_list[k]);

    // timestamp
    out_file_ << measurement_pack_list[k].timestamp_ << "\t"; // pos1 - est

    // output the state vector
    out_file_ << ukf.x_(0) << "\t"; // pos1 - est
    out_file_ << ukf.x_(1) << "\t"; // pos2 - est
    out_file_ << ukf.x_(2) << "\t"; // vel_abs -est
    out_file_ << ukf.x_(3) << "\t"; // yaw_angle -est
    out_file_ << ukf.x_(4) << "\t"; // yaw_rate -est

    // output lidar and radar specific data
    if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER) {
      // sensor type
      out_file_ << "lidar" << "\t";

      // NIS value
      out_file_ << ukf.NIS_laser_ << "\t";

      // output the lidar sensor measurement px and py
      out_file_ << measurement_pack_list[k].raw_measurements_(0) << "\t";
      out_file_ << measurement_pack_list[k].raw_measurements_(1) << "\t";

    } else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR) {
      // sensor type
      out_file_ << "radar" << "\t";

      // NIS value
      out_file_ << ukf.NIS_radar_ << "\t";

      // output radar measurement in cartesian coordinates
      float ro = measurement_pack_list[k].raw_measurements_(0);
      float phi = measurement_pack_list[k].raw_measurements_(1);
      out_file_ << ro * cos(phi) << "\t"; // px measurement
      out_file_ << ro * sin(phi) << "\t"; // py measurement
    }

    // output the ground truth
    out_file_ << gt_pack_list[k].gt_values_(0) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(1) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(2) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(3) << "\n";

    // convert ukf x vector to cartesian to compare to ground truth
    VectorXd ukf_x_cartesian_ = VectorXd(4);

    float x_estimate_ = ukf.x_(0);
    float y_estimate_ = ukf.x_(1);
    float vx_estimate_ = ukf.x_(2) * cos(ukf.x_(3));
    float vy_estimate_ = ukf.x_(2) * sin(ukf.x_(3));

    ukf_x_cartesian_ << x_estimate_, y_estimate_, vx_estimate_, vy_estimate_;

    estimations.push_back(ukf_x_cartesian_);
    ground_truth.push_back(gt_pack_list[k].gt_values_);

  }

  // compute the accuracy (RMSE)
  Tools tools;
  cout << "RMSE" << endl << tools.CalculateRMSE(estimations, ground_truth) << endl;

  // close files
  if (out_file_.is_open()) {
    out_file_.close();
  }

  if (in_file_.is_open()) {
    in_file_.close();
  }

  cout << "Done!" << endl;


  ///// LETS make this an optimization problem with CMAES
  if (false) {
    cmaes_t evo; /* an CMA-ES type struct or "object" */
    double *arFunvals, *const *pop, *xfinal;
    int i;
    int lambda = 7;
    std::vector<double> inxstart = {0.15, 0.5};
    std::vector<double> inrgstddev = {0.1, 0.1};
    VectorXd rubric_gt = VectorXd(4);
    rubric_gt << .07, .08, .35, .25;
    VectorXd ukf_x_cartesian_cmaes = VectorXd(4);


    /* Initialize everything into the struct evo, 0 means default */
    arFunvals = cmaes_init(&evo, (int) inxstart.size(), inxstart.data(), inrgstddev.data(), 0, lambda, "none");
    printf("%s\n", cmaes_SayHello(&evo));
    cmaes_ReadSignals(&evo,
                      "/Users/saminda/Udacity/self_driving_car_term_2/CarND-Unscented-Kalman-Filter-Project/src/cmaes/cmaes_signals.par");  /* write header and initial values */

    /* Iterate until stop criterion holds */
    while (!cmaes_TestForTermination(&evo)) {
      /* generate lambda new search points, sample population */
      pop = cmaes_SamplePopulation(&evo); /* do not change content of pop */

      /* Here we may resample each solution point pop[i] until it
          becomes feasible. function is_feasible(...) needs to be
          user-defined.
          Assumptions: the feasible domain is convex, the optimum is
          not on (or very close to) the domain boundary, initialX is
          feasible and initialStandardDeviations are sufficiently small
          to prevent quasi-infinite looping. */
      /* for (i = 0; i < cmaes_Get(&evo, "popsize"); ++i)
           while (!is_feasible(pop[i]))
             cmaes_ReSampleSingle(&evo, i);
      */

      /* evaluate the new search points using fitfun */
      for (i = 0; i < cmaes_Get(&evo, "lambda"); ++i) {

        //std::cout << "lambda: " << i << std::endl;

        vector<VectorXd> estimations_cmaes;
        vector<VectorXd> ground_truth_cmaes;
        UKF ukf_cmaes;
        double const *pop_x = pop[i];
        ukf_cmaes.std_a_ = std::max(pop_x[0], -pop_x[0]);
        ukf_cmaes.std_yawdd_ = std::max(pop_x[1], -pop_x[1]);

        //std::cout << pop_x[0] << " " << pop_x[1] << std::endl;
        //std:cout << ukf_cmaes.std_a_ << " " << ukf_cmaes.std_yawdd_ << std::endl;

        for (size_t k = 0; k < number_of_measurements; ++k) {
          // Call the UKF-based fusion
          ukf_cmaes.ProcessMeasurement(measurement_pack_list[k]);

          //arFunvals[i] = fitfun(pop[i], (int) cmaes_Get(&evo, "dim"));

          const double xx_estimate_cmaes = ukf_cmaes.x_(0);
          const double yy_estimate_cmaes = ukf_cmaes.x_(1);
          const double vx_estimate_cmaes = ukf_cmaes.x_(2) * cos(ukf_cmaes.x_(3));
          const double vy_estimate_cmaes = ukf_cmaes.x_(2) * sin(ukf_cmaes.x_(3));

          ukf_x_cartesian_cmaes.fill(0.0);
          ukf_x_cartesian_cmaes << xx_estimate_cmaes, yy_estimate_cmaes, vx_estimate_cmaes, vy_estimate_cmaes;

          estimations.push_back(ukf_x_cartesian_cmaes);
          ground_truth.push_back(gt_pack_list[k].gt_values_);

          //std::cout << "\t k: " << k << std::endl;

        }

        Tools tools;
        VectorXd mse = tools.CalculateRMSE(estimations, ground_truth);
        VectorXd res = rubric_gt - mse;
        VectorXd res2 = res.array() * res.array();
        arFunvals[i] = res2.sum();
      }


      /* update the search distribution used for cmaes_SamplePopulation() */
      cmaes_UpdateDistribution(&evo, arFunvals);

      /* read instructions for printing output or changing termination conditions */
      cmaes_ReadSignals(&evo,
                        "/Users/saminda/Udacity/self_driving_car_term_2/CarND-Unscented-Kalman-Filter-Project/src/cmaes/cmaes_signals.par");
      fflush(stdout); /* useful in MinGW */

      xfinal = cmaes_GetNew(&evo, "xbestever"); /* "xbestever" might be used as well */
      std::cout << "xbestever: " << xfinal[0] << " " << xfinal[1] << std::endl;
      free(xfinal);

    }
    printf("Stop:\n%s\n", cmaes_TestForTermination(&evo)); /* print termination reason */
    //cmaes_WriteToFile(&evo, "all", "allcmaes.dat");         /* write final results */

    /* get best estimator for the optimum, xmean */
    xfinal = cmaes_GetNew(&evo, "xmean"); /* "xbestever" might be used as well */

    std::cout << "xmean: " << xfinal[0] << " " << xfinal[1] << std::endl;

    cmaes_exit(&evo); /* release memory */

    /* do something with final solution and finally release memory */
    free(xfinal);
  }

  return 0;
}
