EVENTDIR=../../simu/GRAND/c-aires-reader/src
ANTENNADIR=../../simu/GRAND/c-antenna-test-master/src
libgrand_utillib.a: Makefile aires_util.c hardware_util.c antenna_util.c complex_util.c fft_util.c matrix.c io_util.c reco_util.c
	g++ -c -I$(EVENTDIR) -I$(ANTENNADIR) -I${ROOTSYS}/include/root aires_util.c antenna_util.c hardware_util.c complex_util.c fft_util.c matrix.c reco_util.c io_util.c -std=c++11
	ar -rv libgrand_utillib.a aires_util.o antenna_util.o hardware_util.o complex_util.o matrix.o fft_util.o io_util.o reco_util.o
