import re

with open('src/Libint2IntegralProvider.cpp', 'r') as f:
    text = f.read()

# Remove the old ComputeIntegrals and replace with the three new methods
# We will use regex to find where ComputeIntegrals starts and ends
start_idx = text.find('Integrals Libint2IntegralProvider::ComputeIntegrals() const {')

# The function ends before `IntegralDerivatives Libint2IntegralProvider::ComputeFirstDerivatives() const {`
end_idx = text.find('IntegralDerivatives Libint2IntegralProvider::ComputeFirstDerivatives() const {')

if start_idx != -1 and end_idx != -1:
    new_code = """
T2 Libint2IntegralProvider::ComputeOverlap() const {
    const auto& basis = impl_->basis;
    const auto nbf = static_cast<int>(basis.nbf());
    const auto max_nprim = basis.max_nprim();
    const auto max_l_i = static_cast<int>(basis.max_l());
    libint2::Engine overlap_engine(libint2::Operator::overlap, max_nprim, max_l_i, 0);

    T2 overlap(nbf, nbf);
    overlap.setZero();
    const auto& shell2bf = basis.shell2bf();
    const auto nshell = basis.size();

    for (size_t s1 = 0; s1 < nshell; ++s1) {
        for (size_t s2 = 0; s2 <= s1; ++s2) {
            const auto& shell1 = basis[s1];
            const auto& shell2 = basis[s2];
            const size_t i0 = shell2bf[s1];
            const size_t j0 = shell2bf[s2];
            auto ov_results = overlap_engine.compute(shell1, shell2);
            if (ov_results[0]) {
                FillSymmetricMatrix(overlap, i0, j0, shell1.size(), shell2.size(), ov_results[0]);
            }
        }
    }
    return overlap;
}

T2 Libint2IntegralProvider::ComputeHcore() const {
    const auto& basis = impl_->basis;
    const auto nbf = static_cast<int>(basis.nbf());
    const auto max_nprim = basis.max_nprim();
    const auto max_l_i = static_cast<int>(basis.max_l());
    
    libint2::Engine kinetic_engine(libint2::Operator::kinetic, max_nprim, max_l_i, 0);
    libint2::Engine nuclear_engine(libint2::Operator::nuclear, max_nprim, max_l_i, 0);
    nuclear_engine.set_params(impl_->nuclear_charges);

    T2 hcore(nbf, nbf);
    hcore.setZero();
    const auto& shell2bf = basis.shell2bf();
    const auto nshell = basis.size();

    for (size_t s1 = 0; s1 < nshell; ++s1) {
        for (size_t s2 = 0; s2 <= s1; ++s2) {
            const auto& shell1 = basis[s1];
            const auto& shell2 = basis[s2];
            const size_t i0 = shell2bf[s1];
            const size_t j0 = shell2bf[s2];
            
            auto kin_results = kinetic_engine.compute(shell1, shell2);
            if (kin_results[0]) {
                FillSymmetricMatrix(hcore, i0, j0, shell1.size(), shell2.size(), kin_results[0]);
            }
            auto nuc_results = nuclear_engine.compute(shell1, shell2);
            if (nuc_results[0]) {
                FillSymmetricMatrix(hcore, i0, j0, shell1.size(), shell2.size(), nuc_results[0]);
            }
        }
    }
    return hcore;
}

T4 Libint2IntegralProvider::ComputeERI() const {
    const auto& basis = impl_->basis;
    const auto nbf = static_cast<int>(basis.nbf());
    T4 eri(nbf, nbf, nbf, nbf);
    eri.setZero();

    libint2::Engine coulomb_engine(libint2::Operator::coulomb, basis.max_nprim(), basis.max_l(), 0);
    const auto& shell2bf = basis.shell2bf();
    const auto nshell = basis.size();

    for (size_t s1 = 0; s1 < nshell; ++s1) {
        for (size_t s2 = 0; s2 <= s1; ++s2) {
            for (size_t s3 = 0; s3 <= s1; ++s3) {
                const auto s4_max = (s1 == s3) ? s2 : s3;
                for (size_t s4 = 0; s4 <= s4_max; ++s4) {
                    const auto* buf = coulomb_engine.compute(basis[s1], basis[s2], basis[s3], basis[s4]);
                    if (buf && buf[0] != nullptr) {
                        const size_t i0 = shell2bf[s1];
                        const size_t j0 = shell2bf[s2];
                        const size_t k0 = shell2bf[s3];
                        const size_t l0 = shell2bf[s4];
                        const size_t n1 = basis[s1].size();
                        const size_t n2 = basis[s2].size();
                        const size_t n3 = basis[s3].size();
                        const size_t n4 = basis[s4].size();

                        const double* ptr = buf[0];
                        for (size_t f1 = 0; f1 < n1; ++f1) {
                            const size_t bf1 = i0 + f1;
                            for (size_t f2 = 0; f2 < n2; ++f2) {
                                const size_t bf2 = j0 + f2;
                                for (size_t f3 = 0; f3 < n3; ++f3) {
                                    const size_t bf3 = k0 + f3;
                                    for (size_t f4 = 0; f4 < n4; ++f4) {
                                        const size_t bf4 = l0 + f4;
                                        double val = *ptr++;
                                        eri(bf1, bf2, bf3, bf4) = val;
                                        eri(bf2, bf1, bf3, bf4) = val;
                                        eri(bf1, bf2, bf4, bf3) = val;
                                        eri(bf2, bf1, bf4, bf3) = val;
                                        eri(bf3, bf4, bf1, bf2) = val;
                                        eri(bf3, bf4, bf2, bf1) = val;
                                        eri(bf4, bf3, bf1, bf2) = val;
                                        eri(bf4, bf3, bf2, bf1) = val;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return eri;
}

"""
    text = text[:start_idx] + new_code + text[end_idx:]

with open('src/Libint2IntegralProvider.cpp', 'w') as f:
    f.write(text)

