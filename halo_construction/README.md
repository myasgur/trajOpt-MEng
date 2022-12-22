# Halo Construction
The code here is used to generate perfectly periodic halo orbits given x- and z-amplitudes as well as the parity (Northern/Sourthern). This code is largely based on that included in the appendix of Ethan Geipel's Master's thesis ["HALO Orbit Determination, Integration, and Station Keeping" (2019)](https://docs.wixstatic.com/ugd/7fda7b_1789f8ac40874cac853d3d0a68b5eccc.pdf?index=true). 

I have optimized and extended Ethan's code to generate periodic halo orbits around all of the collinear Lagrange points, not just L1. I also incorporated checks to deal with numerical instabilities, such as inversion of singular matrices and infinite while loops.

## Directory
