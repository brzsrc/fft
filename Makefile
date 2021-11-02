all: fft.ss-test fft.ss-benchmark fft.x86-test fft.x86-benchmark

fft.ss-test: fft.c sincos.c
	/homes/phjk/simplescalar/bin/gcc -o fft.ss-test -DOURSINCOS -DNOFPCONSTANTS -DTEST -DSCALE=3 fft.c sincos.c

fft.x86-test: fft-original.c
	gcc -o fft.x86-test -DTEST fft-original.c -lm

fft.ss-benchmark: fft.c sincos.c
	/homes/phjk/simplescalar/bin/gcc -o fft.ss-benchmark -DOURSINCOS -DNOFPCONSTANTS -DSCALE=12 fft.c sincos.c

fft.x86-benchmark: fft-original.c
	gcc -o fft.x86-benchmark -DSCALE=25 fft-original.c -lm

run-test: fft.ss-test
	/homes/phjk/simplesim/sim-outorder ./fft.ss-test

run-benchmark: fft.ss-benchmark
	/homes/phjk/simplesim/sim-outorder ./fft.ss-benchmark

clean:
	rm -f fft.ss-test fft.ss-benchmark fft.x86-test fft.x86-benchmark fft.o sincos.o 
