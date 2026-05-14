# python_wrapper_thin_sheet_python_wrapper
Running parallel EM simulations for different frequencies in python, modeling conductive sheets and saving results by frequency.

# Python Wrapper

Python wrapper for the Fortran EM thin-sheet modelling code whseet20.par. Used for geophysical simulations of buried conducting sheets.

## How it works

You set up your model in Python, it writes a parameter file, then runs the Fortran binary which does the actual simulation and writes the results.

## Files

- wsheet20.f: main Fortran solver, contains all the EM physics
- SUB_LUD.f: math subroutines needed by wsheet20.f (LU decomposition, numerical integration)
- a.out: the compiled binary, this is what runs.
- wsheet20.par: parameter file that a.out reads, gets generated automatically
- wsheet20.py: the Python wrapper, import from this to set up and run simulations
- python_wrapper_mult_freq.ipynb: main notebook, runs simulations at multiple frequencies in parallel


## Notes

- green_mode=0 computes and saves the LU decomposition, `green_mode=2` reuses a saved one (faster for repeat runs at the same frequency/geometry)
- ColeCole(low, high, tau, alpha) describes the electrical properties. Set `tau=0` for non-polarisable material
- Source types: `plane_wave`, `electric_dipole`, `magnetic_dipole`
