# Fortran compiler
F77 = gfortran-mp-14
# Flags
OPT = -g

# Every *.F and *.f in this directory gets a matching .xOg (basename strips suffix).
SRCS = $(sort $(wildcard *.F) $(wildcard *.f))
EXECUTABLES = $(addsuffix .xOg,$(basename $(SRCS)))

all: $(EXECUTABLES) dex_no_gaia.e

%.xOg: %.F
	$(F77) $(OPT) $< -o $@

%.xOg: %.f
	$(F77) $(OPT) $< -o $@

dex_no_gaia.e: dex_no_gaia.F
	$(F77) $(OPT) $< -o $@

clean:
	rm -f *.xOg dex_no_gaia.e
