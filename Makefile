# it's recommended that you add the following to your startup script:
# export MAKEFLAGS="--check-symlink-times --jobs=$(nproc --all)"
.DEFAULT_GOAL := all

# place in architecture folder, e.g. build/arm64
BUILD_DIR = build/$(shell uname -m)

$(info mpicxx version: $(shell mpicxx --showme:version 2> /dev/null))

ifeq ($(TESTMODE),ON)
else ifeq ($(TESTMODE),OFF)
else ifeq ($(TESTMODE),)
else
  $(error bad value for TESTMODE option. Should be ON or OFF)
endif

ifeq ($(REPRODUCIBLE),ON)
	CXXFLAGS += -DREPRODUCIBLE=true -ffp-contract=off -DEIGEN_DONT_VECTORIZE
	BUILD_DIR := $(BUILD_DIR)_reproducible
	FASTMATH := OFF
else ifeq ($(REPRODUCIBLE),OFF)
else ifeq ($(REPRODUCIBLE),)
else
  $(error bad value for REPRODUCIBLE option. Should be ON or OFF)
endif

CXX := mpicxx
COMPILER_VERSION := $(shell $(CXX) --version)
COMPILER_VERSION_NUMBER := $(shell $(CXX) -dumpversion -dumpfullversion)
COMPILER_VERSION_NUMBER_MAJOR := $(shell echo $(COMPILER_VERSION_NUMBER) | cut -f1 -d.)
$(info $(COMPILER_VERSION))
CXX_STD := c++26

ifneq '' '$(findstring HIP version,$(COMPILER_VERSION))'
	COMPILER_NAME := HIPCC
	CXX_STD := c++23
	CXXFLAGS += -Wno-macro-redefined -Wno-unused-command-line-argument
else ifneq '' '$(findstring clang,$(COMPILER_VERSION))'
	COMPILER_NAME := CLANG
	CXXFLAGS += -Wunsafe-buffer-usage -Wno-unsafe-buffer-usage-in-libc-call -fsafe-buffer-usage-suggestions -Wno-unneeded-internal-declaration
	LDFLAGS += -Wno-unused-command-line-argument

	ifeq '' '$(findstring Apple,$(COMPILER_VERSION))'
		ifeq ($(if $(shell command -v lld),'true','false'), 'true')
			LDFLAGS += -fuse-ld=lld
		endif
	endif
else ifneq '' '$(findstring g++,$(COMPILER_VERSION))'
	COMPILER_NAME := GCC
	# std::stacktrace is available in GCC 14 and later
	# but it is not enabled by default because it slowed down the GitHub CI by > 2x
	ifeq ($(shell expr $(COMPILER_VERSION_NUMBER_MAJOR) \>= 14),1)
		ifeq ($(STACKTRACE),ON)
			CXXFLAGS += -DSTACKTRACE_ON=true -rdynamic
			LDFLAGS += -lstdc++exp
		endif
	endif
	# std=c++26 is not supported on gcc < 14
	ifeq ($(shell expr $(COMPILER_VERSION_NUMBER_MAJOR) \<= 13),1)
		CXX_STD := c++23
	endif
	CXXFLAGS += -Wno-psabi
# 	CXXFLAGS += -Wsuggest-attribute=pure -Wsuggest-attribute=const
else ifneq '' '$(findstring nvc++,$(COMPILER_VERSION))'
	COMPILER_NAME := NVHPC
	CXX_STD := c++23
	# to use the pixi installed libstdc++
# 	CXXFLAGS += --gcc-toolchain=$(PWD)/.pixi/envs/default -Wl,-rpath,$(PWD)/.pixi/envs/default/lib
	ifneq (,$(shell hostname -A | grep .cosma.))
		CXXFLAGS += --gcc-toolchain=/cosma/local/gcc/14.1.0/
	endif

else
	$(warning Unknown compiler)
	COMPILER_NAME := unknown
endif

$(info detected compiler is $(COMPILER_NAME) major version $(COMPILER_VERSION_NUMBER_MAJOR))

CXXFLAGS += -std=$(CXX_STD) -Wall -Wextra -Wpedantic -Wredundant-decls -Wno-unused-parameter -Wsign-compare -Wshadow -isystem third_party

ifneq ($(COMPILER_NAME),NVHPC)
	CXXFLAGS += -Wunused-macros -Werror -Wextra-semi -Wno-unknown-pragmas -Wno-error=cast-function-type -MD -MP -Wno-unused-function
endif

# CXXFLAGS += -DUSE_SIMPSON_INTEGRATOR

# profile-guided optimisation
# generate profile:
# CXXFLAGS += -fprofile-generate="profdataraw"
# for clang, run this to convert the raw data to profdata
# llvm-profdata merge -output=profdata profdataraw/*
# compile with PGO:
# CXXFLAGS += -fprofile-use="profdataraw"

ifeq ($(GPU),ON)
	CXXFLAGS += -DGPU_ON -DUSE_SIMPSON_INTEGRATOR -DBOOST_MATH_NO_EXCEPTIONS -DBOOST_NO_IOSTREAM
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
			CXXFLAGS += -gpu=cc80,rdc
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
			CXXFLAGS += -gpu=cc80,rdc
		else
			CXXFLAGS += -stdpar=multicore
		endif
	else ifeq ($(COMPILER_NAME),HIPCC)
		CXXFLAGS += -fexperimental-library
		ifeq ($(GPU),ON)
			# MI300
			CXXFLAGS += --offload-arch=gfx942 -fgpu-rdc --hipstdpar
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
	# gcc
	# CXXFLAGS += -fopt-info-vec-missed
	# clang
	# CXXFLAGS += -Rpass=loop-vectorize
	# CXXFLAGS += -Rpass-missed=loop-vectorize
	# CXXFLAGS += -Rpass-analysis=loop-vectorize

else
	# sometimes the login nodes have different CPUs
	# to the job nodes. march=native assumes that they are the same
	# (or that the job CPUs support all features of the login CPUs)

	# to get the current CPU architecture, run this:
	# g++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3

	ifneq (,$(shell hostname -A | grep gsi.de))
		# virgo has znver4 nodes in the amd,epyc,9654 feature group
		CXXFLAGS += -march=znver4
		# znver3 in the other partitions and login nodes does not support avx512
		# CXXFLAGS += -march=native -mtune=znver4
	else
		CXXFLAGS += -march=native
	endif

endif

# default to gsl off unless boost or eigen are off
GSL := OFF
ifeq ($(BOOST),OFF)
	CXXFLAGS += -DBOOST_OFF
	GSL := ON
	BUILD_DIR := $(BUILD_DIR)_boostoff
else
	CXXFLAGS += -DBOOST_MATH_STANDALONE
endif

ifeq ($(EIGEN),OFF)
	CXXFLAGS += -DEIGEN_OFF
	GSL := ON
	BUILD_DIR := $(BUILD_DIR)_eigenoff
else
endif

ifeq ($(GSL),ON)
	# GSL (GNU Scientific Library)
	CXXFLAGS += $(shell pkg-config --cflags gsl)
	# Use GSL inline functions
	CXXFLAGS += -DHAVE_INLINE -DGSL_C99_INLINE

	ifeq ($(STATICGSL),)
		# default to static linking
		STATICGSL := ON
	endif

	ifeq ($(STATICGSL),ON)
		gsllibdir := $(shell pkg-config --variable=libdir gsl)
		gsl_objects = $(gsllibdir)/libgsl.a $(gsllibdir)/libgslcblas.a
	else ifeq ($(STATICGSL),OFF)
		LDFLAGS += $(shell pkg-config --libs gsl)
	else
		$(error bad value for STATICGSL option. Should be ON or OFF)
	endif
endif

ifneq ($(MAX_NODE_SIZE),)
	CXXFLAGS += -DMAX_NODE_SIZE=$(MAX_NODE_SIZE)
endif

ifeq ($(TESTMODE),ON)
	CXXFLAGS += -DTESTMODE=true

	CXXFLAGS += -fno-omit-frame-pointer -g

	CXXFLAGS += -D_GLIBCXX_ASSERTIONS
	# CXXFLAGS += -D_GLIBCXX_DEBUG
	# CXXFLAGS += -D_GLIBCXX_DEBUG_BACKTRACE

	# CXXFLAGS += -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_FAST
	CXXFLAGS += -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_EXTENSIVE

	CXXFLAGS += -fsanitize=undefined,address

	BUILD_DIR := $(BUILD_DIR)_testmode
else
	ifeq ($(GSL),ON)
		CXXFLAGS += -DGSL_RANGE_CHECK_OFF
	endif
endif

ifeq ($(OPTIMIZE),OFF)
	BUILD_DIR := $(BUILD_DIR)_optimizeoff
	CXXFLAGS += -O0 -DEIGEN_FAST_MATH=0 -DEIGEN_DONT_ALIGN=1
else
	CXXFLAGS += -O3

	ifeq ($(COMPILER_NAME),CLANG)
		CXXFLAGS += -flto=thin
	else ifeq ($(COMPILER_NAME),GCC)
		CXXFLAGS += -flto=auto -pipe
	endif

	ifeq ($(FASTMATH),OFF)
		BUILD_DIR := $(BUILD_DIR)_nofastmath
		CXXFLAGS += -DEIGEN_FAST_MATH=0 -DEIGEN_DONT_ALIGN=1
	else
		ifeq ($(COMPILER_NAME),NVHPC)
			CXXFLAGS += -fast
		else
			CXXFLAGS += -ffast-math -funsafe-math-optimizations -fno-finite-math-only
		endif
	endif
endif

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
