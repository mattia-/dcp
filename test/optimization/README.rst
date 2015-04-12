This collection of files sets up the case for a test control problem
in which we wish to determine the right inflow condition to reproduce
the flow around a circle with a rectanglle in front of it once the
rectangle is removed

Steps to launch the code:

1) compile:: 

    mkdir build
    cd build
    cmake ..
    make

You may have to adjust the directories where cmake looks for header files and libraries if you installed 
DCP in a non-default location. To do this, use the variables INCLUDE_DIRS and LINK_DIRS when invoking cmake

2) run the first executable

   This executable may take two different command line arguments:
   
   - ``output_file_name`` : 
     the name of the file to which the solution will be saved. 
     It will be a HDF5 file, and format string ``.hdf5`` will be
     automatically appended. 
     The default value is: ``target_solution``
     
   - ``human_readable_print`` : 
     boolean flag. If set to true, also prints velocity and pressure 
     to human readable format
     (in particular, ``pvd`` format - i.e. ``xml`` + ``vtk`` - will be
     used). The name of such files will be the name of the hdf5 output 
     file with the suffix ``_velocity`` or ``_pressure``
     and ``pvd`` extension. The dafault value is ``false``
     
   To launch the executable the command is::
    
    ./src/compute_target_solution --output_file_name <file_name> \ 
    --human_readable_print <boolean_value>

   One or more of the command line options may not be given. In this
   case, defaults values will be used

3) run the second executable

   The command line arguments for this executable are:
   
   - ``target_solution_file_name`` : 
     the name of the file from which the solution will be read. 
     It will be a HDF5 file, and format string ``.hdf5`` will be
     automatically appended. The default value is: ``target_solution``
     
   To launch the executable the command is::
    
    ./src/main --target_solution_file_name <file_name>
