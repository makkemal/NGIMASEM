FORTRAN TEST PROBLEMS
Project Files:
1. usnc.py - Module designed for USNC for curve fitting
2. MainProgram.model - Standard Fortran Code that needs to be written
3. Calibration Environment Jupyter Notebook
4. Test1.inp - Test5.inp
Output File Generated: FOTRANFILE (default=umat.f) and stored in path specified by FORTRANPATH
Software Requirements:
Python Anaconda 3.5: https://www.continuum.io/downloads
-----------------------------------------------------------------------

Two Stages of Fitting of Creep Data:
1. Scalar Material Properties, Youngs Modulus, Linear Dimensional Change and Coefficient of Thermal Expansion, Are Individually Fitted¶
2. Primary and Secondary Creep Models Are Fitted on the Total Creep Strain Data
-----------------------------------------------------------------------
Units are assumed as follows:
1. Temperature degC
2. Fluence (dpa)
3. Force (N)
4. Stress (MPa)
5. Length (mm)
6. Time defined by usnc.GetGetTempPosTimeCode and usnc.GetFluencePosTimeCode as they are the only functions that depend on time
