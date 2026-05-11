FC        = amdflang
OPTIMIZE  ?= yes
DEBUG     ?= no
DEVICE    ?= gpu
ARCH      ?= gfx90a
LAUNCHER  ?=

program = main
source  = main.f90
obj     = $(source:.f90=.o)

FFLAGS := $(EXTRA_FFLAGS) -Wall

ifeq ($(DEBUG),yes)
  FFLAGS += -g
endif

ifeq ($(OPTIMIZE),yes)
  FFLAGS += -O3
endif

ifeq ($(DEVICE),gpu)
  FFLAGS += -fopenmp --offload-arch=$(ARCH)
else
  FFLAGS += -fopenmp
endif

LDFLAGS =

$(program): $(obj)
	$(FC) $(FFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj) output.out

run: $(program)
	$(LAUNCHER) ./$(program) 4 2 2 fixtures/power_4x2 fixtures/temp_4x2 output.out
