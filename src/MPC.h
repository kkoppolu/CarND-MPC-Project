#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 private:
  // predicted x trajectory points
  vector<double> d_trajectory_x;
  // predicted y trajectory points
  vector<double> d_trajectory_y;

  // delta control output from previous solve
  double d_prev_delta;
  // throttle control output from previous solve
  double d_prev_acc;
  // control latency
  double d_latency;

 public:
  MPC(double latency);

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

  // get the trajectory x and y co-ordinates predicted by MPC
  void getTrajectory(vector<double>& trajectoryX, vector<double>& trajectoryY);
};

#endif /* MPC_H */
