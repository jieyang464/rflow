# dmqc — density-matrix quantum chemistry
#
#   make            build the demo driver
#   make test       build and run every unit test
#   make run        build and run the demo driver
#   make clean      remove build artifacts

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
SRCDIR    = src
TESTDIR   = tests
OBJDIR    = build
TARGET    = dmqc

# ── Third-party ──────────────────────────────────────────────────────────────
EIGEN_DIR = third_party/eigen
BASIS_DIR = $(CURDIR)/third_party/libint2_basis

# libint2 is required.  `make` fails with a readable message if it is missing.
LIBINT_CFLAGS := $(shell pkg-config --cflags libint2 2>/dev/null)
LIBINT_LIBS   := $(shell pkg-config --libs   libint2 2>/dev/null)
ifeq ($(strip $(LIBINT_LIBS)),)
  $(error libint2 not found by pkg-config. Install it (Debian/Ubuntu: \
    'sudo apt install libint2-dev') or add its .pc file to PKG_CONFIG_PATH)
endif

# libxc is optional: without it the built-in Slater + VWN5 functionals are used.
LIBXC_CFLAGS := $(shell pkg-config --cflags libxc 2>/dev/null)
LIBXC_LIBS   := $(shell pkg-config --libs   libxc 2>/dev/null)
ifneq ($(strip $(LIBXC_LIBS)),)
  LIBXC_CFLAGS += -DUSE_LIBXC
endif

INCLUDES = -I$(SRCDIR) -I$(EIGEN_DIR) $(LIBINT_CFLAGS) $(LIBXC_CFLAGS)
# Points libint2 at the vendored basis sets so LIBINT_DATA_PATH need not be set
# by hand; an existing LIBINT_DATA_PATH in the environment still wins.
DEFINES  = -DSCFCXX_BASIS_PATH='"$(BASIS_DIR)"'
LDLIBS   = $(LIBINT_LIBS) $(LIBXC_LIBS)

# ── Library sources ──────────────────────────────────────────────────────────
SRCS = \
  $(SRCDIR)/basis.cpp \
  $(SRCDIR)/geometry.cpp \
  $(SRCDIR)/scf.cpp \
  $(SRCDIR)/fock_builders.cpp \
  $(SRCDIR)/jk_builder.cpp \
  $(SRCDIR)/IDensityUpdater.cpp \
  $(SRCDIR)/SzaboHeHIntegral.cpp \
  $(SRCDIR)/Libint2IntegralProvider.cpp \
  $(SRCDIR)/Libint2DerivativeProvider.cpp \
  $(SRCDIR)/effective_densities.cpp \
  $(SRCDIR)/gradient.cpp \
  $(SRCDIR)/energy_gradient.cpp \
  $(SRCDIR)/geometry_optimizer.cpp \
  $(SRCDIR)/grid/molecular_grid.cpp \
  $(SRCDIR)/dft/dft_helper.cpp \
  $(SRCDIR)/xc/lsda.cpp \
  $(SRCDIR)/xc/libxc_wrapper.cpp \
  $(SRCDIR)/ci.cpp

OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))

TESTS = \
  test_integral_providers \
  test_density_updater \
  test_dft \
  test_gradient \
  test_geometry_optimizer

TEST_BINS = $(addprefix $(OBJDIR)/, $(TESTS))

# ── Rules ────────────────────────────────────────────────────────────────────
.PHONY: all test run clean check-basis

all: $(TARGET)

$(TARGET): $(OBJDIR)/main.o $(OBJS)
	$(CXX) $(CXXFLAGS) $^ $(LDLIBS) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c $< -o $@

$(OBJDIR)/%: $(TESTDIR)/%.cpp $(OBJS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -I$(TESTDIR) $(DEFINES) $< $(OBJS) $(LDLIBS) -o $@

# Each suite runs to completion even if an earlier one fails, so a single run
# reports every failure rather than only the first.
test: check-basis $(TEST_BINS)
	@fail=0; \
	for binary in $(TEST_BINS); do \
	  ./$$binary || fail=1; \
	done; \
	if [ $$fail -eq 0 ]; then echo "ALL SUITES PASSED"; else echo "SOME SUITES FAILED"; fi; \
	exit $$fail

run: $(TARGET)
	./$(TARGET)

check-basis:
	@if [ ! -f "$(BASIS_DIR)/sto-3g.g94" ]; then \
	  echo "Basis sets not found in $(BASIS_DIR)."; \
	  echo "Run ./install_basis_sets.sh first."; \
	  exit 1; \
	fi

clean:
	rm -rf $(OBJDIR) $(TARGET)
