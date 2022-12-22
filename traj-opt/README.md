# Trajectory Optimization
This code: 
- sets up the three-body system (Sun, Earth, spacecraft) and relevant parameters (spacecraft wet mass, thrust, Isp, time of flight, canonical units, etc.)
- defines the IVP: search for the costates (Lagrange multipliers) that 
    1. satisfy the terminal constraints of the originial TPBVP (refer to [Woollands & Eggl (2020)](https://trs.jpl.nasa.gov/bitstream/handle/2014/52355/CL%2320-0708.pdf?sequence=1&isAllowed=y)), and 
    2. minimize the cost function (maximize the spacecraft mass at the end of the orbital transfer)
- sets up the full variational equations of motion (equations of motion with optimal control + state transition matrix dynamics)
- *attempts* to optimize the performance of the indirect single shooting procedure by
    1. determining the minimum timestep for the fixed-step integrator `ode4` to maintain equivalent precision to higher-order adaptive time-step integrators such as `ode78`
    2. constraining the search space of initial costates in R7 using a uniformly-distributed random process and feasibility tolerance to hone in on regions where viable solutions may exist (currently flawed; in development)

<!-- ## Directory
- `main.m`  -->
