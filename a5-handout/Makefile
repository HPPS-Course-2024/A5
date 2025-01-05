CC?=cc
CFLAGS?=-Wall -Wextra -pedantic -lm -fopenmp -O3
DEBUG?=-g
OPENMP?=-fopenmp
PROGRAMS= genparticles particles2text warnings2text cmpparticles nbody nbody-bh

all: $(PROGRAMS)

clean:
	rm -f $(PROGRAMS) *.o

%: %.o util.o
	$(CC) -o $@ $*.o util.o $(CFLAGS)

%.o: %.c
	$(CC) -c $< $(CFLAGS) $(DEBUG) $(OPENMP)

.SUFFIXES: # Disables implicit rules.  (Don't worry about it.)
