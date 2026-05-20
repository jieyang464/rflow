CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra
XC_CXXFLAGS ?=
XC_LDLIBS ?=

# ─── Directory layout ────────────────────────────────────────────────────────

SRCDIR   = src
OBJDIR   = build
TARGET   = test_bin

# ─── Third-party paths ──────────────────────────────────────────────────────

# Eigen (header-only, vendored in third_party/)
EIGEN_DIR = third_party/eigen

# Basis-set data 
BASIS_DIR  = third_party/libint2_basis

# ─── Aggregate compiler / linker flags ───────────────────────────────────────

FEATURE_DEFS = -DSCFCXX_ENABLE_LIBINT2_PROVIDER=0
INCLUDES = -I$(SRCDIR) -I$(EIGEN_DIR)
LDFLAGS  =
LDLIBS   = $(XC_LDLIBS)

# ─── Source & object lists ───────────────────────────────────────────────────

SRCS = $(SRCDIR)/test.cpp \
       $(SRCDIR)/scf.cpp \
       $(SRCDIR)/fock_builders.cpp \
       $(SRCDIR)/jk_builder.cpp \
       $(SRCDIR)/ci.cpp \
       $(SRCDIR)/SzaboHeHIntegral.cpp \
       $(SRCDIR)/IDensityUpdater.cpp \
       $(SRCDIR)/xc/vxc_evaluator.cpp \
       $(SRCDIR)/xc/libxc_wrapper.cpp \
       $(SRCDIR)/xc/vxc_libxc_grid.cpp

OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))

# =============================================================================
#  Primary targets
# =============================================================================

# Default: just compile the project.
all: $(TARGET)

# Link objects into the final binary.
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(FEATURE_DEFS) $(OBJS) $(LDFLAGS) $(LDLIBS) -o $@

# Compile each .cpp → .o, creating build/ if necessary.
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(FEATURE_DEFS) $(XC_CXXFLAGS) $(INCLUDES) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

# Build then run the binary.
run: $(TARGET)
	./$(TARGET)

# =============================================================================
#  Cleanup
# =============================================================================

clean:
	 rm -rf $(OBJDIR) $(TARGET)

clean-all: clean
	 rm -rf $(BASIS_DIR)

.PHONY: all clean clean-all run
