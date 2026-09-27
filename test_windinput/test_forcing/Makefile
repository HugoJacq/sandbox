CFLAGS += -I$(DOCUMENT_ROOT)/sandbox -O2
CFLAGS += -disable-dimensions # skip dimensions check for now

all: no_forcing.tst no_forcing/plots

linear_wave_wind_input.tst: CFLAGS+= -DFORCING=1 -DHOLD_FORCING=1

no_forcing.c: linear_wave_wind_input.c
	ln -sf linear_wave_wind_input.c no_forcing.c

exact_forcing.c: linear_wave_wind_input.c
	ln -sf linear_wave_wind_input.c exact_forcing.c

exact_forcing.tst: CFLAGS += -DEXACT_FORCING=1




include $(BASILISK)/Makefile.defs


