# Gamma #
### Generic Synthesis C++ Library


# About #

Gamma is a cross-platform, C++ library for doing generic synthesis and 
filtering of signals. It contains helpful mathematical functions, 
types, such as vectors and complex numbers, an assortment of sequence 
generators, and many other objects for signal processing tasks. 
It is oriented towards real-time sound and graphics synthesis, but is 
equally useful for non-real-time tasks.


# Compilation Instructions #

The source code can either be built into a library or directly compiled from source into an application. In the following, the base directory is where this `README` file is located.

## Building a Library

### Make (Linux, OS X, mingw)
In most cases, simply running

	make

will build the library with automatically detected platform settings. See `Makefile.config` for other build options.

There are several other rules within Makefile. These are:

	make			- builds static library
	make install		- installs library into DESTDIR
	make clean		- removes binaries from build folder
	make test		- performs unit tests

The script `run.sh` can be used to compile and run examples and other source files against the Gamma library. For example,

	./run.sh examples/oscillator/sine.cpp

To only compile the source file without running, include AUTORUN=0 after the source file. Binaries are located in the automatically generated `build/` directory. On OSX, the Gamma library will be linked to the pre-compiled dependent libraries in `external/lib_osx`. On Linux, use `apt-get` to install the necessary dependent libraries.


### Xcode (OS X)
You may also build the library using the supplied Xcode project.
1. Open `project/xcode/gamma.xcodeproj`
2. Build the target `libgamma{.a, .dylib}`. The library will be in project build folder.


## Compiling Directly From Source
Gamma can easily be compiled directly from source into an existing project.

Make sure to pass in the following flags to the compiler:

	-D__STDC_CONSTANT_MACROS
	-finline-functions (or -O3)
	-fpeel-loops

## Dependencies

PortAudio is required ONLY if you are using Gamma's AudioIO class (defined in Gamma/AudioIO.h). If you do not wish to use audio i/o, then pass the flag

	NO_AUDIO_IO=1

into make or, if not using make, exclude `src/AudioIO.cpp` from your project.

libsndfile may be used as the backend for the SoundFile class by passing the flag

	USE_LIBSNDFILE=1

into make.

## License

Gamma is distributed under a permissive free software license. Please see the LICENSE file for details.