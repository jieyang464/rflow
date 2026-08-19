CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

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

INCLUDES = -I$(SRCDIR) -I$(EIGEN_DIR)
LDFLAGS  =
LDLIBS   =

# ─── Source & object lists ───────────────────────────────────────────────────

SRCS = $(SRCDIR)/test.cpp \
       $(SRCDIR)/scf.cpp \
       $(SRCDIR)/fock_builders.cpp \
       $(SRCDIR)/jk_builder.cpp \
       $(SRCDIR)/ci.cpp \
       $(SRCDIR)/SzaboHeHIntegral.cpp \
       $(SRCDIR)/IDensityUpdater.cpp

# To enable Libint2 support manually, add the provider source and library:
# SRCS += $(SRCDIR)/Libint2IntegralProvider.cpp
# LDLIBS += -lint2

OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))

# =============================================================================
#  Primary targets
# =============================================================================

# Default: just compile the project.
all: $(TARGET)

# Link objects into the final binary.
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LDLIBS) -o $@

# Compile each .cpp → .o, creating build/ if necessary.
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

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
