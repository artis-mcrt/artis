# it's recommended that you add the following to your startup script:
# export MAKEFLAGS="--check-symlink-times --jobs=$(nproc --all)"
.DEFAULT_GOAL := all

# place in architecture folder, e.g. build/arm64
BUILD_DIR = build/$(shell uname -m)

CXX := mpicxx
$(info mpicxx version: $(shell mpicxx --showme:version 2> /dev/null))

ifeq ($(TESTMODE),ON)
else ifeq ($(TESTMODE),OFF)
else ifeq ($(TESTMODE),)
else
  $(error bad value for TESTMODE option. Should be ON or OFF)
endif

ifeq ($(REPRODUCIBLE),ON)
	CXXFLAGS += -DREPRODUCIBLE=true -ffp-contract=off
	BUILD_DIR := $(BUILD_DIR)_reproducible
else ifeq ($(REPRODUCIBLE),OFF)
else ifeq ($(REPRODUCIBLE),)
else
  $(error bad value for REPRODUCIBLE option. Should be ON or OFF)
endif

COMPILER_VERSION := $(shell $(CXX) --version)
$(info $(COMPILER_VERSION))
ifneq '' '$(findstring clang,$(COMPILER_VERSION))'
  COMPILER_NAME := CLANG
  CXXFLAGS += -flto=thin
else ifneq '' '$(findstring g++,$(COMPILER_VERSION))'
  COMPILER_NAME := GCC
  CXXFLAGS += -flto=auto
else ifneq '' '$(findstring nvc++,$(COMPILER_VERSION))'
  COMPILER_NAME := NVHPC
else
  $(warning Unknown compiler)
  COMPILER_NAME := unknown
endif

$(info detected compiler is $(COMPILER_NAME))

ifeq ($(COMPILER_NAME),NVHPC)
	CXXFLAGS += -std=c++20
else
	CXXFLAGS += -std=c++23 -ftree-vectorize -Wunused-macros -Werror -Wno-error=unknown-pragmas -MD -MP -ftrivial-auto-var-init=pattern
endif

CXXFLAGS += -fstrict-aliasing -DMPI_ON=true
# CXXFLAGS += -DUSE_SIMPSON_INTEGRATOR=true

# profile-guided optimisation
# generate profile:
# CXXFLAGS += -fprofile-generate="profdataraw"
# for clang, run this to convert the raw data to profdata
# llvm-profdata merge -output=profdata profdataraw/*
# compile with PGO:
# CXXFLAGS += -fprofile-use="profdataraw"

ifeq ($(GPU),ON)
	CXXFLAGS += -DGPU_ON=true -DUSE_SIMPSON_INTEGRATOR=true
	BUILD_DIR := $(BUILD_DIR)_gpu
else ifeq ($(GPU),OFF)
else ifeq ($(GPU),)
else
    $(error bad value for GPU option. Should be ON or OFF)
endif

ifeq ($(OPENMP),ON)
  ifeq ($(STDPAR),ON)
    $(error cannot combine OPENMP and STDPAR)
  endif
  BUILD_DIR := $(BUILD_DIR)_openmp

	ifeq ($(COMPILER_NAME),NVHPC)
		ifeq ($(GPU),ON)
			CXXFLAGS += -mp=gpu -gpu=mem:unified
		else
			CXXFLAGS += -mp
		endif
	else ifeq ($(COMPILER_NAME),CLANG)
		CXXFLAGS += -Xpreprocessor -fopenmp
		LDFLAGS += -lomp
	else ifeq ($(COMPILER_NAME),GCC)
		CXXFLAGS += -fopenmp
	endif

else ifeq ($(OPENMP),OFF)
else ifeq ($(OPENMP),)
else
    $(error bad value for OPENMP option. Should be ON or OFF)
endif

ifeq ($(STDPAR),ON)
  CXXFLAGS += -DSTDPAR_ON=true
  BUILD_DIR := $(BUILD_DIR)_stdpar

  ifeq ($(COMPILER_NAME),NVHPC)
		ifeq ($(GPU),ON)
			CXXFLAGS += -stdpar=gpu -gpu=mem:unified
			# CXXFLAGS += -gpu=cc80
		else
			CXXFLAGS += -stdpar=multicore
		endif
  else ifeq ($(COMPILER_NAME),CLANG)
		CXXFLAGS += -fexperimental-library
		# LDFLAGS += -ltbb
		# LDFLAGS += -Xlinker -debug_snapshot
  else ifeq ($(COMPILER_NAME),GCC)
		LDFLAGS += -ltbb
  endif

else ifeq ($(STDPAR),OFF)
else ifeq ($(STDPAR),)
else
  $(error bad value for STDPAR option. Should be ON or OFF)
endif

ifneq ($(STDPAR),ON)
	# triggers errors for the onedpl headers
	CXXFLAGS += -Wundef
endif

ifeq ($(shell uname -s),Darwin)
# 	macOS

	ifeq ($(shell uname -m),arm64)
#	 	On Arm, -mcpu combines -march and -mtune
		CXXFLAGS += -mcpu=native
	else
#		On x86, -march implies -mtune
		CXXFLAGS += -march=native
	endif

	CXXFLAGS += -fno-omit-frame-pointer -g
#	CXXFLAGS += -Rpass=loop-vectorize
#	CXXFLAGS += -Rpass-missed=loop-vectorize
#	CXXFLAGS += -Rpass-analysis=loop-vectorize

	# add -lprofiler for gperftools
	# LDFLAGS += $(LIB)
	# LDFLAGS += -lprofiler

else
	# sometimes the login nodes have slightly different CPUs
	# to the job nodes. Try to find the lowest common denominator here
	# to enable vector extensions
	# CXXFLAGS += -march=cascadelake
	# CXXFLAGS += -march=skylake-avx512
	# CXXFLAGS += -march=icelake-server

	# to get the current CPU architecture, run this:
	# g++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3
	ifneq (,$(shell hostname -A | grep gsi.de))
		# virgo has some AMD nodes and some Intel.
		# As of 6th July 2024, the login nodes are cascadelake, and we select zen3 nodes to run on.
		CXXFLAGS += -march=x86-64-v3 -mtune=znver3
	else
		# for GitHub actions, checksums must match with different assigned CPUs, so avoid -march=native (use lowest common denominator)
		# update: all zenver3 now?

		CXXFLAGS += -march=native
	endif

endif

# GSL (GNU Scientific Library)
CXXFLAGS += $(shell pkg-config --cflags gsl)

ifeq ($(STATICGSL),)
	# default to dynamic linking
	STATICGSL := OFF
endif

ifeq ($(STATICGSL),ON)
	gsllibdir := $(shell pkg-config --variable=libdir gsl)
	gsl_objects = $(gsllibdir)/libgsl.a $(gsllibdir)/libgslcblas.a
	BUILD_DIR := $(BUILD_DIR)_staticgsl
else ifeq ($(STATICGSL),OFF)
	LDFLAGS += $(shell pkg-config --libs gsl)
else
    $(error bad value for STATICGSL option. Should be ON or OFF)
endif

# Use GSL inline functions
CXXFLAGS += -DHAVE_INLINE -DGSL_C99_INLINE

ifneq ($(MAX_NODE_SIZE),)
	CXXFLAGS += -DMAX_NODE_SIZE=$(MAX_NODE_SIZE)
endif

ifeq ($(TESTMODE),ON)
	CXXFLAGS += -DTESTMODE=true -D_LIBCPP_DEBUG=0

	CXXFLAGS += -D_GLIBCXX_ASSERTIONS
	# CXXFLAGS += -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_BACKTRACE=1

	CXXFLAGS +=  -fno-omit-frame-pointer -g

	# CXXFLAGS += -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_EXTENSIVE
	CXXFLAGS += -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_DEBUG

	CXXFLAGS += -fsanitize=undefined,address

	BUILD_DIR := $(BUILD_DIR)_testmode
else
	# skip GSL range checking for better performance
	CXXFLAGS += -DTESTMODE=false -DGSL_RANGE_CHECK_OFF
endif

ifeq ($(OPTIMIZE),OFF)
	BUILD_DIR := $(BUILD_DIR)_optimizeoff
	CXXFLAGS += -O0
else
	CXXFLAGS += -O3
	ifeq ($(FASTMATH),OFF)
		BUILD_DIR := $(BUILD_DIR)_nofastmath
		CXXFLAGS += -Wno-unknown-pragmas
	else
		ifeq ($(COMPILER_NAME),NVHPC)
			CXXFLAGS += -fast
		else
			CXXFLAGS += -ffast-math -funsafe-math-optimizations -fno-finite-math-only -fopenmp-simd
		endif
	endif
endif

CXXFLAGS += -Winline -Wall -Wpedantic -Wredundant-decls -Wno-unused-parameter -Wno-unused-function -Wno-inline -Wsign-compare -Wshadow

# sn3d.cc and exspec.cc have main() defined
common_files := $(filter-out sn3d.cc exspec.cc, $(wildcard *.cc))

sn3d_files = $(common_files) sn3d.cc
sn3d_objects = $(addprefix $(BUILD_DIR)/,$(sn3d_files:.cc=.o))
sn3d_dep = $(sn3d_objects:%.o=%.d)

exspec_files = $(common_files) exspec.cc
exspec_objects = $(addprefix $(BUILD_DIR)/,$(exspec_files:.cc=.o))
exspec_dep = $(exspec_objects:%.o=%.d)

.ONESHELL:
define version_h
constexpr const char* GIT_VERSION = \"$(shell git describe --dirty --always --tags)\";
constexpr const char* GIT_BRANCH = \"$(shell git symbolic-ref --short HEAD 2>/dev/null || git rev-parse --short HEAD )\";
constexpr const char* GIT_STATUS = \"$(shell git status --short)\";
endef

$(shell echo "$(version_h)" > version_tmp.h)
$(shell test -f version.h || touch version.h)

ifneq ($(shell cat version.h),$(shell cat version_tmp.h))
  $(info updating version.h)
  $(shell mv version_tmp.h version.h)
else
  $(shell rm version_tmp.h)
endif

$(shell mkdir -p $(BUILD_DIR))

$(shell echo "$(COMPILER_VERSION)" > $(BUILD_DIR)/compiler_tmp.txt)
$(shell test -f $(BUILD_DIR)/compiler.txt || touch $(BUILD_DIR)/compiler.txt)
ifneq ($(shell cat $(BUILD_DIR)/compiler.txt),$(shell cat $(BUILD_DIR)/compiler_tmp.txt))
  $(info detected compiler change)
  $(shell mv $(BUILD_DIR)/compiler_tmp.txt $(BUILD_DIR)/compiler.txt)
else
  $(shell rm $(BUILD_DIR)/compiler_tmp.txt)
endif

all: sn3d exspec

$(BUILD_DIR)/%.o: %.cc Makefile $(BUILD_DIR)/compiler.txt
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

check: $(sn3d_files)
	run-clang-tidy $(sn3d_files)

$(BUILD_DIR)/sn3d: $(sn3d_objects)
	$(CXX) $(CXXFLAGS) $(sn3d_objects) $(gsl_objects) $(LDFLAGS) -o $(BUILD_DIR)/sn3d
-include $(sn3d_dep)

sn3d: $(BUILD_DIR)/sn3d
	ln -sf $(BUILD_DIR)/sn3d sn3d

$(BUILD_DIR)/sn3dwhole: $(sn3d_files) version.h artisoptions.h Makefile $(BUILD_DIR)/compiler.txt
	$(CXX) $(CXXFLAGS) $(sn3d_files) $(gsl_objects) $(LDFLAGS) -o $(BUILD_DIR)/sn3dwhole
-include $(sn3d_dep)

sn3dwhole: $(BUILD_DIR)/sn3dwhole
	ln -sf $(BUILD_DIR)/sn3dwhole sn3d

$(BUILD_DIR)/exspec: $(exspec_objects)
	$(CXX) $(CXXFLAGS) $(exspec_objects) $(gsl_objects) $(LDFLAGS) -o $(BUILD_DIR)/exspec
-include $(exspec_dep)

exspec: $(BUILD_DIR)/exspec
	ln -sf $(BUILD_DIR)/exspec exspec

.PHONY: clean sn3d sn3dwhole exspec

clean:
	rm -rf sn3d exspec build *.o *.d
