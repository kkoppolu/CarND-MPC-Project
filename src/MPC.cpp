#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

namespace {
static const size_t N = 12;

// fg = [cost, x, y, v, psi, cte, epsi, delta, acc] - each having N elements
// vars = [x, y, v, psi, cte, epsi, delta, acc] - each having N elements
static const int x_start = 0;
static const int y_start = x_start + N;
static const int psi_start = y_start + N;
static const int v_start = psi_start + N;
static const int cte_start = v_start + N;
static const int epsi_start = cte_start + N;
static const int delta_start = epsi_start + N;
static const int a_start = delta_start + N - 1;
}  // namespace

class FG_eval {
 private:
  // reference velocity
  static const double ref_v;

 public:
  // timestep delta
  static const double dt;

  // This value assumes the model presented in the classroom is used.
  //
  // It was obtained by measuring the radius formed by running the vehicle in
  // the simulator around in a circle with a constant steering angle and
  // velocity on a flat terrain.
  //
  // Lf was tuned until the the radius formed by the simulating the model
  // presented in the classroom matched the previous radius.
  //
  // This is the length from front to CoG that has a similar radius.
  static const double Lf;

  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable
    // values (state & actuators) NOTE: You'll probably go back and forth
    // between this function and the Solver function below.
    fg[0] = 0;

    // minimize CTE and epsi
    for (size_t t = 0; t < N; t++) {
      fg[0] += 30 * CppAD::pow(vars[cte_start + t], 2);
      fg[0] += 15 * CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2);

      // Minimize the use of actuators.
      for (size_t t = 0; t < N - 1; t++) {
        fg[0] += CppAD::pow(vars[delta_start + t], 2);
        fg[0] += CppAD::pow(vars[a_start + t], 2);
      }

      // Minimize the value gap between sequential actuations.
      for (size_t t = 0; t < N - 2; t++) {
        fg[0] +=
            400 *
            CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
        fg[0] += 100 * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
      }

      // initialization - offset by 1 for cost function
      fg[1 + x_start] = vars[x_start];
      fg[1 + y_start] = vars[y_start];
      fg[1 + psi_start] = vars[psi_start];
      fg[1 + v_start] = vars[v_start];
      fg[1 + cte_start] = vars[cte_start];
      fg[1 + epsi_start] = vars[epsi_start];
    }

    // model predictions
    for (size_t t = 1; t < N; t++) {
      // The state at time t+1 .
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];

      // The state at time t.
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];

      AD<double> f0 = coeffs[0] + coeffs[1] * x0 +
                      coeffs[2] * CppAD::pow(x0, 2) +
                      coeffs[3] * CppAD::pow(x0, 3);
      AD<double> f0_dot = coeffs[1] + (2 * coeffs[2] * x0) +
                          (3 * coeffs[3] * CppAD::pow(x0, 2));
      AD<double> psides0 = CppAD::atan(f0_dot);

      // Recall the equations for the model:
      // x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
      // y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
      // psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
      // v_[t] = v[t-1] + a[t-1] * dt
      // cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
      // epsi[t] = psi[t-1] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
      fg[1 + cte_start + t] =
          cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + epsi_start + t] =
          epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }
  }
};

const double FG_eval::dt = 0.05;
const double FG_eval::Lf = 2.67;
const double FG_eval::ref_v = 50;

//
// MPC class definition implementation.
//
MPC::MPC(double latency) : d_prev_delta(0), d_prev_acc(0), d_latency(latency) {}

MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  size_t n_vars = 6 * N + 2 * (N - 1);
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (size_t i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // Set the initial variable values taking into account the controls latency
  double f0 = coeffs[0] + coeffs[1] * x + coeffs[2] * CppAD::pow(x, 2) +
              coeffs[3] * CppAD::pow(x, 3);
  double f0_dot =
      coeffs[1] + (2 * coeffs[2] * x) + (3 * coeffs[3] * CppAD::pow(x, 2));
  double psides = CppAD::atan(f0_dot);
  const double delt = d_latency / 1000;
  vars[x_start] = x + v * CppAD::cos(psi) * delt;
  vars[y_start] = y + v * CppAD::sin(psi) * delt;
  vars[psi_start] = psi + v * d_prev_delta / FG_eval::Lf * delt;
  vars[v_start] = v + d_prev_acc * delt;
  vars[cte_start] = (f0 - y) + (v * CppAD::sin(epsi) * delt);
  vars[epsi_start] =
      (psi - psides) + v * d_prev_delta / FG_eval::Lf * delt;  // epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (size_t i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (size_t i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (size_t i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (size_t i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  // auto cost = solution.obj_value;
  // std::cout << "Cost " << cost << std::endl;

  d_prev_delta = solution.x[delta_start];
  d_prev_acc = solution.x[a_start];
  d_trajectory_x.clear();
  d_trajectory_y.clear();
  for (size_t i = 0; i < N; ++i) {
    size_t xIdx = x_start + i;
    size_t yIdx = y_start + i;
    d_trajectory_x.push_back(solution.x[xIdx]);
    d_trajectory_y.push_back(solution.x[yIdx]);
  }
  // creates a 2 element double vector.
  return {d_prev_delta, d_prev_acc};
}

void MPC::getTrajectory(vector<double>& trajectoryX,
                        vector<double>& trajectoryY) {
  trajectoryX.clear();
  trajectoryY.clear();

  std::copy(d_trajectory_x.begin(), d_trajectory_x.end(),
            back_inserter(trajectoryX));
  std::copy(d_trajectory_y.begin(), d_trajectory_y.end(),
            back_inserter(trajectoryY));
}