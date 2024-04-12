# Photo-Voltaic CuRve AnalyZEr (PV-CRAZE)
The PV-CRAZE library is developed to work with descriptive current/voltage (IV) parameters for solar cells and modules. In particular, PV-CRAZE determines the extended solar cell parameters (ESP) [1]. To this end it implements robust methods to analyze experimental IV measurements, and determin ethe ESP. In addition, PV-CRAZE can fit various curve models to the ESP to reconstruct an approximate IV from its ESP. Supported models are among others:

- One Diode model
- Karmalka-Haneefa model[2] 
- PV-Spline [1]
- PV-Bézier [1]

The code also includes methods to use the ESP in conjunction with the IEC60891 norm to transpose the ESP as a function of irradiance and temperature. Experimental code is provided to fit the IEC60891 model to a set of ESP at different irradiances and temperatures (your milage may vary).

Bindings to C, python, and GNU Octave are provided.

## Installing
To install use the usual

    ./configure
    make
	sudo make install
	
This installs the pvcraze library on your system. You may need to run ldconfig so programs can find it:

	sudo ldconfig

If that does not work, run it with the library install directory, e.g. 

	sudo ldconfig /usr/local/lib
	
To make it permanent you may need to add the library install directory to /etc/ld.so.conf.d/ 

If you want to install the octave interface use

    octave --eval "pkg install -verbose src/octave-pv-craze.tar.gz"

For the python interaface:

    pip install src/

You can test the library and the octave and python interfaces by running the tests.

    cd tests
    ./runtest.sh

## Usage

See DOC.md for details on the library. You can also look at the tests in the tests
directory for examples.

 
## References
[1] B.Pieters, "Extended Solar Cell Parameters - General Purpose Descriptive I/V Parameters for Solar Cells." to be submitted (2023)

[2] S. Karmalkar and S. Haneefa, "A physically based explicit J-V model of a solar cell for simple design calculations." IEEE Electron Device Letters 29, no. 5 (2008): 449-451.
 
