This program (both Fortran version and Python version) simulates a Browning motion.
There should be a proportional relationship between MSD(dt) and dt.

There are 4 arguments you can assign, use argument "--help" only to see details.

The straight-forward way of calculating MSD is a O(N^2) task which may be really 
time-consuming, especially for the Python version. Hence, a fft (fast Fourier transform) 
method is used, which reduced the time complexity to O(NlogN).

For the Fortran version, the program is linked to a fftw library (http://fftw.org/).
This program is tested under Windows environment, and the windows version of fftw 
can be downloaded from http://fftw.org/install/windows.html. For Linux users just 
download the source code from the official website and compile fftw yourself. Note 
that the "Makefile" file need to be slightly modified to link fftw to this program
correctly.

The standard normal distribution random numbers in the Fortran version is generated
by an advanced version of Box-Muller's algorithm.

