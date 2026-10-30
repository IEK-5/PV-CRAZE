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

The Python package is built independently of the autotools build: it
compiles the library sources into its extension module, so it does not need
the installed library. It does need a C compiler and the GSL and LAPACK
development files. From the repository root:

    pip install .            # or: pip install .[pandas]

You can test the library and the octave and python interfaces by running the tests.

    cd tests
    ./runtest.sh

The Python tests alone run with `python3 -m pytest tests/python`.

## Usage

Python:

```python
import pvcraze as pc

res = pc.analyze(V, I, models="PV-Bezier")   # V, I: 1-D arrays of one IV curve
res.ok, res.flag().messages                  # analysis sanity checks
res.esp.voc, res.esp.isc, res.esp.ff         # extended solar cell parameters
model = res.fits["PV-Bezier"].params
model.iv([0.0, 0.3])                         # evaluate the curve model
```

Batches work the same way: pass a list of `(n, 2)` arrays and every field
becomes an array. Everything has docstrings (`help(pc.analyze)`); DOC.md
has a Python section with an overview and a migration table from the old
dictionary based API.

See DOC.md for details on the library and the Octave interface. You can also
look at the tests in the tests directory for examples.

 
## References
[1] B.Pieters, "Extended Solar Cell Parameters - General Purpose Descriptive I/V Parameters for Solar Cells." to be submitted (2023)

[2] S. Karmalkar and S. Haneefa, "A physically based explicit J-V model of a solar cell for simple design calculations." IEEE Electron Device Letters 29, no. 5 (2008): 449-451.
 
