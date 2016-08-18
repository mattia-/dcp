# check on CMake version
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

ENABLE_TESTING ()

ADD_SUBDIRECTORY (algebraic_problems)
ADD_SUBDIRECTORY (expressions_and_subdomains)
ADD_SUBDIRECTORY (linear_differential_problems)
ADD_SUBDIRECTORY (nonlinear_differential_problems)
ADD_SUBDIRECTORY (optimization)
ADD_SUBDIRECTORY (time_dependent_problems)

ADD_TEST (NAME algebraic_problems COMMAND algebraic_problems)
ADD_TEST (NAME expressions_and_subdomains COMMAND expressions_and_subdomains)
ADD_TEST (NAME linear_problems COMMAND linear_problem)
ADD_TEST (NAME nonlinear_problems COMMAND nonlinear_problem)
ADD_TEST (NAME optimization_external_load COMMAND optimization_external_load_control)
ADD_TEST (NAME time_dependent_convection_diffusion COMMAND time_dependent_convection_diffusion)
ADD_TEST (NAME time_dependent_cavity_NS COMMAND time_dependent_cavity_NS)
