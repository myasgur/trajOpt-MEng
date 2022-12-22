This code is heavily based on that included in the appendix of Ethan Geipel's Master's thesis "HALO Orbit Determination, Integration, and Station Keeping" (2019). 
Here is a link to the pdf of said thesis: https://docs.wixstatic.com/ugd/7fda7b_1789f8ac40874cac853d3d0a68b5eccc.pdf?index=true

I have optimized Ethan's code to generate periodic halo orbits around all of the collinear Lagrange points, not just L1. I also incorporated checks to deal with numerical instabilities, such as inversion of singular matrices or infinite while loops.
