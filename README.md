# scfcpp

A minimal Hartree–Fock (RHF) playground using Eigen for linear algebra and Libint2 for Gaussian integrals. Includes routines to compute overlap, kinetic, nuclear attraction, and two-electron integrals (ERIs), plus a simple SCF harness.

## Requirements

- A C++17 compiler (g++/clang++)
- pkg-config
- Eigen 3.4+
- Libint2 2.7+

### Ubuntu/Debian (apt)

```bash
sudo apt update
sudo apt install -y build-essential pkg-config libeigen3-dev libint2-dev
```

If `libint2-dev` isn’t found, enable the universe repo and retry:
```bash
sudo add-apt-repository universe
sudo apt update
sudo apt install -y libint2-dev
```

### Fedora (dnf)

```bash
sudo dnf install -y eigen3-devel libint-devel pkgconf-pkg-config
```

### Arch (pacman)

```bash
sudo pacman -S --needed eigen libint2 pkgconf
```

### Build from source (fallback)

- Eigen (header-only): https://gitlab.com/libeigen/eigen
- Libint2 releases: https://github.com/evaleev/libint/releases

Example (installs to `$HOME/.local`):
```bash
# Eigen
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
tar xzf eigen-3.4.0.tar.gz && cd eigen-3.4.0
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/.local -DBUILD_TESTING=OFF
cmake --build build --target install

# Libint2
cd ..
wget https://github.com/evaleev/libint/releases/download/v2.7.2/libint-2.7.2.tgz
tar xzf libint-2.7.2.tgz && cd libint-2.7.2
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/.local -DBUILD_SHARED_LIBS=ON -DINSTALL_GTEST=OFF
cmake --build build -j
cmake --install build

# If installed to a non-system prefix, you might need this at link or runtime:
export PKG_CONFIG_PATH=$HOME/.local/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$HOME/.local/lib:$LD_LIBRARY_PATH
```

Verify installs:
```bash
pkg-config --modversion eigen3
pkg-config --modversion libint2
```

## Build

Use the provided Makefile (uses pkg-config):
```bash
make
```
This builds `build/scfcpp_demo` by compiling sources in `src/` and tests in `test/`.

Alternatively, compile manually:
```bash
g++ -std=c++17 -O2 \
  $(pkg-config --cflags eigen3 libint2) \
  src/*.cpp test/*.cpp \
  -o build/scfcpp_demo \
  $(pkg-config --libs libint2) -pthread
```

If pkg-config isn’t available, try:
```bash
g++ -std=c++17 -O2 \
  -I/usr/include/eigen3 \
  src/*.cpp test/*.cpp \
  -o build/scfcpp_demo \
  -lint2 -pthread
```

## Run

```bash
./build/scfcpp_demo
```
This will:
- Print the number of basis functions (nbf)
- Print the core Hamiltonian (Hcore)
- Print non-zero ERI elements for a small demo molecule (HeH+ with sto-3g)

## Basis sets

The project includes many `.g94` basis files under `basis/`. The code attempts to use `./basis/<name>.g94` automatically by setting `LIBINT_DATA_PATH=./basis` internally when a matching file exists. You can also pass an explicit `basis_file_path` to `Molecule`.

## Layout

- `src/` — library code (HF and Molecule)
- `test/` — small demo/test drivers
  - `eri_hcore_demo.cpp` prints Hcore and ERIs for a toy molecule
- `basis/` — Gaussian basis set files
- `build/` — build outputs (created by Makefile)

## Troubleshooting

- Linker can’t find libint2 at runtime:
  - Set runtime path on link: `-Wl,-rpath,$HOME/.local/lib`
  - Or export: `export LD_LIBRARY_PATH=$HOME/.local/lib:$LD_LIBRARY_PATH`
- `pkg-config` can’t find libint2:
  - Ensure `PKG_CONFIG_PATH` includes your install prefix: `export PKG_CONFIG_PATH=$HOME/.local/lib/pkgconfig:$PKG_CONFIG_PATH`
- `fatal error: Eigen/Dense: No such file or directory`:
  - Install `libeigen3-dev` (Debian/Ubuntu) or add `-I/path/to/eigen3`.

## License

TBD