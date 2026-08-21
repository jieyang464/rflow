# scfcxx

A self-directed implementation of **density-matrix-based SCF**: an unrestricted
Hartree–Fock and Kohn–Sham program in which the density matrix is converged by
*rotating it* rather than by diagonalizing a Fock matrix.

This is a **proof of concept**. The claim it makes good on is narrow and
specific: *swapping the Fock diagonalization for a commutator-driven density
update converges to the same density matrix.* It is not a fast program, and
[Scope and limitations](#scope-and-limitations) says plainly what is missing.

---

## The idea

The usual SCF cycle builds a Fock matrix and solves the generalized eigenvalue
problem `FC = SCε` to get a new density. That diagonalization is `O(N³)` and
cannot be made to scale linearly, so it eventually becomes the wall that
linear-scaling methods run into.

The alternative is to treat the density matrix itself as the variable, and to
move it along the manifold of valid densities with an exponential
parameterization. Following Helgaker, Jørgensen and Olsen (Ch. 10, §10.7), the
density is updated by a similarity transform in the overlap metric:

```
D  ←  e^(−XS)  D  e^(SX)
```

with the generator taken from the **S-metric commutator** of the Fock and
density matrices:

```
X  =  h · (F D S  −  S D F)
```

Three properties make this work.

**It is the steepest-descent direction.** `R = FDS − SDF` is the gradient of the
energy with respect to an orbital rotation, and it vanishes exactly at a
stationary point — which is why it is also the DIIS error matrix. Because `F`,
`D` and `S` are symmetric, `R` is *antisymmetric*, so `Tr(R·R) = −‖R‖²_F < 0`
and a positive step `h` lowers the energy. (The sign is easy to get backwards:
with the opposite sign the iteration converges just as happily, to the *highest*
stationary point.)

**It preserves the constraints exactly.** `e^(−XS) D e^(SX)` is a similarity
transform in the `S` metric, so both the electron count `Tr(DS)` and the
idempotency `DSD = D` are conserved — no orthogonalization, no re-imposition of
occupation numbers.

**It needs no fock matrix diagonalization.** The transform is evaluated by a
Baker–Campbell–Hausdorff expansion in the S-metric adjoint action,

```
e^(−XS) D e^(SX)  =  Σₖ (−1)ᵏ/k!  adᵏ(D),      ad(A) = X S A − A S X
```

so each iteration costs only matrix multiplications. Truncating the series
breaks idempotency slightly, which is repaired by one McWeeny purification:

```
D  ←  3 D S D  −  2 D S D S D
```

Nothing in the iteration diagonalizes anything. A Fock matrix is diagonalized
exactly once, *after* convergence, and only to produce canonical orbitals for
whatever comes next — a deliberate design point, controlled by
`SCFSettings::finalize_orbitals`.

## Main archetecture features
The project applies abstract base class interfaces extensively on places like
integral (derivative) providers, density matrice updator, Fock builder effective 
density builders to achieve a modular and maintainable design.
Such as, it currently provides two integral providers one based on the Szabo 
book's HeH+ numerical example, another being the libint2, and is flexible with 
adding new integral providers
---

## Installing

**Requirements**

| | |
|---|---|
| C++17 compiler | g++ or clang++ |
| [libint2](https://github.com/evaleev/libint) | required — molecular integrals |
| Eigen 3 | vendored in `third_party/eigen`, nothing to install |
| [libxc](https://libxc.gitlab.io/) | optional — built-in Slater + VWN5 is used without it |

**Debian / Ubuntu**

```bash
sudo apt install build-essential pkg-config libint2-dev
sudo apt install libxc-dev          # optional
```

**Basis sets.** libint2 reads basis sets from `.g94` files at runtime. Fetch
them into the repository once:

```bash
./install_basis_sets.sh
```

This writes to `third_party/libint2_basis/`, and the build hard-codes that path
so `LIBINT_DATA_PATH` does not need to be set by hand. If you do export
`LIBINT_DATA_PATH`, it takes precedence.

**Build and run**

```bash
make            # build the demonstration driver
make test       # build and run all five test suites
make run        # converge HeH+ and H2O both ways, side by side
```

`make test` runs every suite to completion even if one fails, so a single run
reports all failures rather than only the first.


## Scope and limitations

This is a proof of concept, and the following are deliberate:

**No integral screening.** This is the big one. Genuine linear scaling requires
Schwarz screening to skip negligible shell quartets, and none is applied. The
density update is `O(N³)` matrix multiplication and the J/K build is `O(N⁴)`, so
the program demonstrates that the *algorithm* works without demonstrating that
it *scales*. `ForEachEriDerivative` is shaped to make screening a natural
addition — it is a decision about a shell quartet, and that is the level the
provider iterates at.


**Initial guess on the density is still done by Hcore diagonalization.** To 
achieve full linear scaling, the initial guess must be done without matrix 
daigonalization, but that is not implemented yet.


### Planned

- **Schwarz screening** 
- **Initial guess without diagonalization**
- **Sparse density matrices from Eigen**
- **MP2 with on-the-fly t2** 
- **MP2 gradient entirely based on AO integrals**
- **CCD and CIS(D_infinite)**


---

## References

The method:

1. T. Helgaker, P. Jørgensen and J. Olsen, *Molecular Electronic-Structure
   Theory*, Wiley (2000), Ch. 10, §10.7 — exponential parameterization and
   direct optimization of the density matrix.
2. R. McWeeny, *Some recent advances in density matrix theory*,
   Rev. Mod. Phys. **32**, 335 (1960) — the purification transform.

Reference data and components:

3. A. Szabo and N. S. Ostlund, *Modern Quantum Chemistry*, Dover (1996),
   §3.5.2 — the HeH⁺ integrals used as the fixed reference.
4. E. F. Valeev, *Libint: a library for the evaluation of molecular integrals of
   many-body operators over Gaussian functions*, https://libint.valeyev.net/
5. A. D. Becke, *A multicenter numerical integration scheme for polyatomic
   molecules*, J. Chem. Phys. **88**, 2547 (1988) — the DFT grid.
6. P. A. M. Dirac, Proc. Cambridge Philos. Soc. **26**, 376 (1930) — local
   exchange.
7. S. H. Vosko, L. Wilk and M. Nusair, Can. J. Phys. **58**, 1200 (1980) —
   the VWN correlation fit (parameterization V).
8. P. Pulay, *Ab initio calculation of force constants and equilibrium
   geometries*, Mol. Phys. **17**, 197 (1969) — the basis-set-derivative term
   in the gradient.

---

## Licence

MIT. See [LICENSE](LICENSE).
