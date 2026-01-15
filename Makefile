CXX ?= g++
CCACHE ?=

# Build mode: release (default) or debug
MODE ?= release
ifeq ($(MODE),debug)
  CXXOPT := -O0 -g
else
  CXXOPT := -O2 -DNDEBUG
endif

INCLUDES := -Isrc $(shell pkg-config --cflags eigen3 libint2)
CXXFLAGS ?= -std=c++17 $(CXXOPT) $(INCLUDES) -MMD -MP
LDFLAGS ?= $(shell pkg-config --libs libint2) -pthread

# Optional sparse matrix build: pass SPARSE=1 to use Eigen::SparseMatrix
ifeq ($(SPARSE),1)
	CXXFLAGS += -DSCFCPP_USE_EIGEN_SPARSE
endif

# Optional LibXC support: pass XC=1 to enable and link libxc
ifeq ($(XC),1)
	CXXFLAGS += -DUSE_LIBXC
	LDFLAGS += $(shell pkg-config --libs libxc)
endif

# Core sources (exclude legacy/experimental translation units)
# - Exclude src/main.cpp (standalone demo)
# - Exclude legacy HF/KS and experimental commutator_scf sources
SRCS := $(filter-out src/main.cpp src/hf.cpp src/ks.cpp src/commutator_scf.cpp,$(wildcard src/*.cpp))

# Object and dependency files
OBJDIR := build/obj
OBJS := $(patsubst src/%.cpp,$(OBJDIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)

# Static library for core code to enable incremental builds
LIB := build/libscfcpp.a

TESTS := build/eri_hcore_demo build/scf_compare build/scf_mixed_demo build/trace_init build/uhf_density_log build/dft_hf_equiv build/szabo_heh_example build/commutator_bch_demo

.PHONY: all clean build run

all: $(TESTS)

build:
	mkdir -p build $(OBJDIR)

# Compile sources to objects with dependency tracking
$(OBJDIR)/%.o: src/%.cpp | build
	$(CCACHE) $(CXX) $(CXXFLAGS) -c $< -o $@

# Archive static library from object files
$(LIB): $(OBJS) | build
	@echo "AR $@"
	ar rcs $@ $(OBJS)

# Link tests against the static library (fast relink)
build/eri_hcore_demo: $(LIB) test/eri_hcore_demo.cpp | build
	$(CCACHE) $(CXX) $(CXXFLAGS) test/eri_hcore_demo.cpp $(LIB) -o $@ $(LDFLAGS)

build/scf_compare: $(LIB) test/scf_compare.cpp | build
	$(CCACHE) $(CXX) $(CXXFLAGS) test/scf_compare.cpp $(LIB) -o $@ $(LDFLAGS)

build/scf_mixed_demo: $(LIB) test/scf_mixed_demo.cpp | build
	$(CCACHE) $(CXX) $(CXXFLAGS) test/scf_mixed_demo.cpp $(LIB) -o $@ $(LDFLAGS)

build/trace_init: $(LIB) test/trace_init.cpp | build
	$(CCACHE) $(CXX) $(CXXFLAGS) test/trace_init.cpp $(LIB) -o $@ $(LDFLAGS)

build/uhf_density_log: $(LIB) test/uhf_density_log.cpp | build
	$(CCACHE) $(CXX) $(CXXFLAGS) test/uhf_density_log.cpp $(LIB) -o $@ $(LDFLAGS)

build/dft_hf_equiv: $(LIB) test/dft_hf_equiv.cpp | build
	$(CCACHE) $(CXX) $(CXXFLAGS) test/dft_hf_equiv.cpp $(LIB) -o $@ $(LDFLAGS)

build/szabo_heh_example: $(LIB) test/szabo_heh_example.cpp | build
	$(CCACHE) $(CXX) $(CXXFLAGS) test/szabo_heh_example.cpp $(LIB) -o $@ $(LDFLAGS)

build/commutator_bch_demo: $(LIB) test/commutator_bch_demo.cpp | build
	$(CCACHE) $(CXX) $(CXXFLAGS) test/commutator_bch_demo.cpp $(LIB) -o $@ $(LDFLAGS)

run: all
	./build/scf_compare

clean:
	rm -rf build

# Include auto-generated dependency files
-include $(DEPS)
