# Low-Thrust Trajectory Optimization 
This repository contains the code for my Master of Engineering project: "Design of Fuel-Optimal Low-Thrust Trajectories to Service the James Webb Space Telescope" (2022). 

## Directory
- `halo_construction` contains the necessary files to generate a nominal periodic halo orbit used for setting terminal conditions on the trajectory optimization. An initial third order approximation is first generated using [Richardson's (1980)](https://link.springer.com/content/pdf/10.1007/BF01229511.pdf?pdf=button) approach which is then iteratively adjusted using [Howell's (1984)](https://www.proquest.com/docview/303183159?pq-origsite=gscholar&fromopenview=true) differential correction method to ensure periodicity.

- `traj-opt` contains the main optimization code which sets up the corresponding initial value problem (IVP) for the two-point boundary value problem (TPBVP) defined by [Woollands & Eggl (2020)](https://trs.jpl.nasa.gov/bitstream/handle/2014/52355/CL%2320-0708.pdf?sequence=1&isAllowed=y).

- `util` contains miscellaneous utility functions to offload and organize the primary project code
