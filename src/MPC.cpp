#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

/************* Step 4: Set and experiment timesteps and durations  ***************** 
* Choose N = 10 dt = 0.1 (s)
* Other N values has been tested. 
* N = 5 will not work.
* N = 20 will not work either.
* N = 15 works as well as N = 10. 
* dt = 0.1s is chosen, because it is easier to deal with the 100ms latency. 
************************************************************************************/

//Set the timestep length and duration
size_t N = 10;
double dt = 0.1;

//Set the targets for cte, epsi, and speed 
double ref_cte = 0;
double ref_epsi = 0;
double ref_v = 50 * 0.44704; //convert speed from MPH to m/s;

//Set the indices for variable and constraints
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;


class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    /************* Step 5: Set Fg vector ************************************************ 
    * fg[0] is cost
    * For each state variable, the value for the starting timepoint is the initial value.
    * such as fg[1 + x_start] = vars[x_start];
    * The values for the rest of N-1 timepoints are constraints based on motion model.
    * such as  fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
    ************************************************************************************/
   
    // COST
    // Note the values before each equations are weights.
    fg[0] = 0;
    for ( int t = 0; t <N; t++) {
       fg[0] += 1000*CppAD::pow(vars[cte_start + t] - ref_cte, 2);
       fg[0] += 1000*CppAD::pow(vars[epsi_start + t] - ref_epsi, 2);
       fg[0] += 1 * CppAD::pow(vars[v_start + t] - ref_v, 2);
    }

    for ( int t = 0; t < N-1; t++ ) {
        fg[0] += 100 * CppAD::pow(vars[delta_start + t], 2);
        fg[0] += 10 * CppAD::pow(vars[a_start + t], 2);
    }

    for (int t = 0; t < N-2; t++) {
        fg[0] += 2000 * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
        fg[0] += 10 *CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }

    // Initial state vars: 6
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // Set up the constrains based on motion model 
    for ( int t = 1; t < N; t++) {

        // states at time t 
        AD<double> x0 = vars[x_start + t - 1];
        AD<double> y0 = vars[y_start + t - 1];
        AD<double> psi0 = vars[psi_start + t - 1];
        AD<double> v0 = vars[v_start + t - 1];
        AD<double> cte0 = vars[cte_start + t - 1];
        AD<double> epsi0 = vars[epsi_start + t - 1];

        // states at time t + 1
        AD<double> x1 = vars[x_start + t];
        AD<double> y1 = vars[y_start + t];
        AD<double> psi1 = vars[psi_start + t];
        AD<double> v1 = vars[v_start + t];
        AD<double> cte1 = vars[cte_start + t];
        AD<double> epsi1 = vars[epsi_start + t];

        // controls, we just need N-1 controls for N time points.
        AD<double> delta0 = vars[delta_start + t - 1];
        AD<double> a0 = vars[a_start + t - 1];

        // targets for y and psi, calculated by using the passed coeffs, which are fitted by way points
        AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;; 
        AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3*coeffs[3]*x0*x0); 

        // constraints.
        // Note that the signs for delta0 need extra attentions.
        fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
        fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
        fg[1 + psi_start + t] = psi1 - (psi0 - v0 * delta0 / Lf * dt);
        fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
        fg[1 + cte_start + t] = cte1 - ( (y0 - f0) - v0 * CppAD::sin(epsi0) * dt);
        fg[1 + epsi_start + t] = epsi1 - ( (psi0 - psides0) - v0 * delta0 / Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  // size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  /************* Step 6: Set the model variables and constraints ********************* 
  * The state is a 6 element vector: x, y, psi, v, cte, epsi 
  * note that steer_latency and throttle_latency are just for the first 100 ms.
  * N timesteps: 6*N + 2* (N-1) variables
  ************************************************************************************/
  double x = state[0];  
  double y = state[1];  
  double psi = state[2];  
  double v = state[3];  
  double cte = state[4];  
  double epsi = state[5];  

  double steer_latency = state[6];  
  double throttle_latency = state[7];  
  
  size_t n_vars = N*6 + (N-1)*2;
  size_t n_constraints = N*6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);

  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0.0;
  }

  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  // use current steer and throttle to cover the latency period.
  vars[delta_start] = steer_latency;
  vars[a_start] = throttle_latency;

  //Set lower and upper limits for variables.
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // Set initial rate of change for psi, from current steer angle
  // For only one time point, because the setting for dt is 100ms. and the latency is also 100ms
  vars_lowerbound[delta_start] = steer_latency;
  vars_upperbound[delta_start] = steer_latency;
  // Set rate of change for psi for all other timepoints
  for (int i = delta_start + 1; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332 * Lf;
    vars_upperbound[i] = 0.436332 * Lf;
  }

  // Set initial throttle, from current throttle 
  // For only one time point, because the setting for dt is 100ms. and the latency is also 100ms
  vars_lowerbound[a_start] = throttle_latency;
  vars_upperbound[a_start] = throttle_latency;
  // Set throttle for all other timepoints
  for (int i = a_start + 1; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }


  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);

  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // constraints for the initial states
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
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  vector<double> result;

  // Note that we return the [delta_start + 1] and [a_start + 1] values. 
  // because the 100ms is just 1 interval for the timepoints.
  // we can just return the second timepoint values as actuator commands.
  // these commands will spend 100ms to take action.
  result.push_back(solution.x[delta_start + 1]);
  result.push_back(solution.x[a_start + 1]);

  // for display the green line.
  for (int i = 0; i < N -1; i++) {
    result.push_back(solution.x[x_start + i + 1]);
    result.push_back(solution.x[y_start + i + 1]);
  }
  
  return result; 
}
