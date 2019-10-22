/******************************************************************************
Copyright 2019 The Regents of the University of California.
All Rights Reserved.

Permission to copy, modify and distribute any part of this Software for
educational, research and non-profit purposes, without fee, and without
a written agreement is hereby granted, provided that the above copyright
notice, this paragraph and the following three paragraphs appear in all
copies.

Those desiring to incorporate this Software into commercial products or
use for commercial purposes should contact the:
Office of Innovation & Commercialization
University of California, San Diego
9500 Gilman Drive, Mail Code 0910
La Jolla, CA 92093-0910
Ph: (858) 534-5815
FAX: (858) 534-7345
E-MAIL: invent@ucsd.edu

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR
DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY
OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS. THE UNIVERSITY OF CALIFORNIA MAKES NO
REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR
EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE
SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
******************************************************************************/

#include "potential/dispersion/dispersion.h"
#include "potential/electrostatics/helpme.h"

namespace disp {

void Dispersion::Initialize(const std::vector<double> sys_c6_long_range, const std::vector<double> &sys_xyz,
                            const std::vector<std::string> &mon_id, const std::vector<size_t> &num_atoms,
                            const std::vector<std::pair<std::string, size_t> > &mon_type_count,
                            const bool do_grads = true, const std::vector<double> &box = {}, std::vector<double> *virial) {
    sys_c6_long_range_ = sys_c6_long_range;
    sys_xyz_ = sys_xyz;
    mon_id_ = mon_id;
    num_atoms_ = num_atoms;
    mon_type_count_ = mon_type_count;
    do_grads_ = do_grads;
    box_ = box;
    box_inverse_ = box.size() ? InvertUnitCell(box) : std::vector<double>{};
    use_pbc_ = box.size();

    natoms_ = sys_c6_long_range_.size();
    size_t natoms3 = 3 * natoms_;
    phi_ = std::vector<double>(natoms_, 0.0);
    xyz_ = std::vector<double>(natoms3, 0.0);
    grad_ = std::vector<double>(natoms3, 0.0);
    virial_ = std::vector<double>(9,0.0);
    sys_grad_ = std::vector<double>(natoms3, 0.0);
    c6_long_range_ = std::vector<double>(natoms_, 0.0);
    sys_phi_ = std::vector<double>(natoms_, 0.0);
    
    calc_virial_ = true;

    if (virial == 0) { // set virial_ to 0 in case we dont want to calculate it
        calc_virial_ = false;
    }

    ReorderData();
}

void Dispersion::SetNewParameters(const std::vector<double> &xyz, bool do_grads = true, const double cutoff = 100.0,
                                  const std::vector<double> &box = {}) {
    sys_xyz_ = xyz;
    box_ = box;
    box_inverse_ = box.size() ? InvertUnitCell(box) : std::vector<double>{};
    use_pbc_ = box.size();
    do_grads_ = do_grads;
    cutoff_ = cutoff;
    std::fill(grad_.begin(), grad_.end(), 0.0);
    std::fill(phi_.begin(), phi_.end(), 0.0);
    std::fill(virial_.begin(), virial_.end(),0.0);

    ReorderData();
}

void Dispersion::ReorderData() {
    // Organize xyz so we have
    // x1_1 x1_2 ... y1_1 y1_2... z1_1 z1_2 ... x2_1 x2_2 ...
    // where xN_M is read as coordinate x of site N of monomer M
    // for the first monomer type. Then follows the second, and so on.
    // As an example, for a water dimer (OHHM), the system xyz would be
    // sys_xyz = {Ox, Oy, Oz, Hx, Hy, Hz, Hx, Hy, ...}
    // And after reordering,
    // xyz = {Ox, Ox, Oy, Oy, Oz, Oz, Hx, Hx, Hy, Hy, Hz, Hz, Hx, Hx, ... ...}

    size_t fi_mon = 0;
    size_t fi_crd = 0;
    size_t fi_sites = 0;
    for (size_t mt = 0; mt < mon_type_count_.size(); mt++) {
        size_t ns = num_atoms_[fi_mon];
        size_t nmon = mon_type_count_[mt].second;
        size_t nmon2 = nmon * 2;
        for (size_t m = 0; m < nmon; m++) {
            size_t mns = m * ns;
            size_t mns3 = mns * 3;
            for (size_t i = 0; i < ns; i++) {
                size_t inmon = i * nmon;
                size_t inmon3 = 3 * inmon;
                xyz_[inmon3 + m + fi_crd] = sys_xyz_[fi_crd + mns3 + 3 * i];
                xyz_[inmon3 + m + fi_crd + nmon] = sys_xyz_[fi_crd + mns3 + 3 * i + 1];
                xyz_[inmon3 + m + fi_crd + nmon2] = sys_xyz_[fi_crd + mns3 + 3 * i + 2];
                c6_long_range_[fi_sites + m + inmon] = sys_c6_long_range_[fi_sites + mns + i];
            }
        }
        fi_mon += nmon;
        fi_sites += nmon * ns;
        fi_crd += nmon * ns * 3;
    }
}

double Dispersion::GetDispersion(std::vector<double> &grad,std::vector<double> *virial) {
    CalculateDispersion();

    size_t fi_mon = 0;
    size_t fi_crd = 0;
    size_t fi_sites = 0;

    if (virial != 0) {
        (*virial)[0] += virial_[0];
        (*virial)[1] += virial_[1];
        (*virial)[2] += virial_[2];
        (*virial)[3] += virial_[3];
        (*virial)[4] += virial_[4];
        (*virial)[5] += virial_[5];
        (*virial)[6] += virial_[6];
        (*virial)[7] += virial_[7];
        (*virial)[8] += virial_[8];
    }

    for (size_t mt = 0; mt < mon_type_count_.size(); mt++) {
        size_t ns = num_atoms_[fi_mon];
        size_t nmon = mon_type_count_[mt].second;
        size_t nmon2 = nmon * 2;
        for (size_t m = 0; m < nmon; m++) {
            size_t mns = m * ns;
            size_t mns3 = mns * 3;
            for (size_t i = 0; i < ns; i++) {
                size_t inmon = i * nmon;
                size_t inmon3 = 3 * inmon;

                sys_phi_[fi_sites + mns + i] = phi_[fi_sites + m + inmon];

                grad[fi_crd + mns3 + 3 * i] += grad_[inmon3 + m + fi_crd];
                grad[fi_crd + mns3 + 3 * i + 1] += grad_[inmon3 + m + fi_crd + nmon];
                grad[fi_crd + mns3 + 3 * i + 2] += grad_[inmon3 + m + fi_crd + nmon2];
            }
        }
        fi_mon += nmon;
        fi_sites += nmon * ns;
        fi_crd += nmon * ns * 3;
    }

    return disp_energy_;
}

void Dispersion::CalculateDispersion() {
    disp_energy_ = 0.0;
    std::fill(phi_.begin(), phi_.end(), 0.0);
    // Max number of monomers
    size_t maxnmon = mon_type_count_.back().second;
    // Parallelization
    size_t nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel  // omp_get_num_threads() needs to be inside
                      // parallel region to get number of threads
    {
        if (omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
    }
#endif

    // This part looks at sites inside the same monomer
    // Reset first indexes
    size_t fi_mon = 0;
    size_t fi_crd = 0;
    size_t fi_sites = 0;

    // Excluded sets
    excluded_set_type exc12;
    excluded_set_type exc13;
    excluded_set_type exc14;

    // Loop over each monomer type
    for (size_t mt = 0; mt < mon_type_count_.size(); mt++) {
        size_t ns = num_atoms_[fi_mon];
        size_t nmon = mon_type_count_[mt].second;
        size_t nmon2 = 2 * nmon;

        // Obtain excluded pairs for monomer type mt
        systools::GetExcluded(mon_id_[fi_mon], exc12, exc13, exc14);

        // For parallel region
        std::vector<std::vector<double> > phi_pool;
        std::vector<std::vector<double> > grad_pool;
        std::vector<double> energy_pool(nthreads,0.0);
        std::vector<std::vector<double> > virial_pool;


        for (size_t i = 0; i < nthreads; i++) {
            phi_pool.push_back(std::vector<double>(nmon * ns,0.0));
            grad_pool.push_back(std::vector<double>(nmon * ns * 3,0.0));
            virial_pool.push_back(std::vector<double>(9,0.0));
        }
        // Loop over each pair of sites
        for (size_t i = 0; i < ns - 1; i++) {
            size_t inmon = i * nmon;
            size_t inmon3 = inmon * 3;
            for (size_t j = i + 1; j < ns; j++) {
                // Continue only if i and j are not bonded
                bool is12 = systools::IsExcluded(exc12, i, j);
                bool is13 = systools::IsExcluded(exc13, i, j);
                bool is14 = systools::IsExcluded(exc14, i, j);
                double disp_scale_factor = (is12 || is13 || is14) ? 0 : 1;
                double c6, d6;
                double c6i = c6_long_range_[fi_sites + i * nmon];
                double c6j = c6_long_range_[fi_sites + j * nmon];
                GetC6(mon_id_[fi_mon], mon_id_[fi_mon], i, j, c6, d6);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
                for (size_t m = 0; m < nmon; m++) {
                    int rank = 0;
#ifdef _OPENMP
                    rank = omp_get_thread_num();
#endif
                    double p1[3], g1[3];
                    double phi_i = 0.0;
                    p1[0] = xyz_[fi_crd + inmon3 + m];
                    p1[1] = xyz_[fi_crd + inmon3 + nmon + m];
                    p1[2] = xyz_[fi_crd + inmon3 + nmon2 + m];
                    std::fill(g1, g1 + 3, 0.0);
                    energy_pool[rank] += disp6(c6, d6, c6i, c6j, p1, xyz_.data() + fi_crd, g1, grad_pool[rank].data(), phi_i,
                                          phi_pool[rank].data(), nmon, nmon, m, m + 1, i, j, disp_scale_factor,
                                          do_grads_, cutoff_, ewald_alpha_, box_, box_inverse_, &virial_pool[rank]);

                    grad_pool[rank][inmon3 + m] += g1[0];
                    grad_pool[rank][inmon3 + nmon + m] += g1[1];
                    grad_pool[rank][inmon3 + nmon2 + m] += g1[2];
                    phi_pool[rank][inmon + m] += phi_i;
                }
            }
        }

        // Compress data in phi and grad and virial
        for (size_t rank = 0; rank < nthreads; rank++) {
            size_t kend = phi_pool[rank].size();
            for (size_t k = 0; k < kend; k++) {
                phi_[fi_sites + k] += phi_pool[rank][k];
            }
            kend = grad_pool[rank].size();
            for (size_t k = 0; k < kend; k++) {
                grad_[fi_crd + k] += grad_pool[rank][k];
            }
            for (size_t k = 0; k < 9; k++) {
                virial_[k] += virial_pool[rank][k];
            }
            disp_energy_ += energy_pool[rank];
        }

        fi_mon += nmon;
        fi_sites += nmon * ns;
        fi_crd += nmon * ns * 3;
    }

    // Sites corresponding to different monomers
    // Declaring first indexes
    size_t fi_mon1 = 0;
    size_t fi_sites1 = 0;
    size_t fi_mon2 = 0;
    size_t fi_sites2 = 0;
    size_t fi_crd1 = 0;
    size_t fi_crd2 = 0;

    // Loop over all monomer types
    for (size_t mt1 = 0; mt1 < mon_type_count_.size(); mt1++) {
        size_t ns1 = num_atoms_[fi_mon1];
        size_t nmon1 = mon_type_count_[mt1].second;
        size_t nmon12 = nmon1 * 2;
        fi_mon2 = fi_mon1;
        fi_sites2 = fi_sites1;
        fi_crd2 = fi_crd1;

        // For each monomer type mt1, loop over all the other monomer types
        // mt2 >= mt1 to avoid double counting
        for (size_t mt2 = mt1; mt2 < mon_type_count_.size(); mt2++) {
            size_t ns2 = num_atoms_[fi_mon2];
            size_t nmon2 = mon_type_count_[mt2].second;

            // Check if monomer types 1 and 2 are the same
            // If so, same monomer won't be done, since it has been done in
            // previous loop.
            bool same = (mt1 == mt2);

            // For parallel region
            std::vector<std::vector<double> > phi1_pool;
            std::vector<std::vector<double> > phi2_pool;
            std::vector<std::vector<double> > grad1_pool;
            std::vector<std::vector<double> > grad2_pool;
            std::vector<double> energy_pool(nthreads,0.0);
            std::vector<std::vector<double> > virial_pool;
            for (size_t i = 0; i < nthreads; i++) {
                phi1_pool.push_back(std::vector<double>(nmon1 * ns1,0.0));
                phi2_pool.push_back(std::vector<double>(nmon2 * ns2,0.0));
                grad1_pool.push_back(std::vector<double>(nmon1 * ns1 * 3,0.0));
                grad2_pool.push_back(std::vector<double>(nmon2 * ns2 * 3,0.0));
                virial_pool.push_back(std::vector<double>(9,0.0));
            }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for (size_t m1 = 0; m1 < nmon1; m1++) {
                int rank = 0;
#ifdef _OPENMP
                rank = omp_get_thread_num();
#endif
                size_t m2init = same ? m1 + 1 : 0;
                for (size_t i = 0; i < ns1; i++) {
                    size_t inmon1 = i * nmon1;
                    size_t inmon13 = inmon1 * 3;
                    double c6i = c6_long_range_[fi_sites1 + i * nmon1];
                    std::vector<double> xyz_sitei(3);
                    std::vector<double> g1(3, 0.0);
                    double phi_i = 0.0;
                    xyz_sitei[0] = xyz_[fi_crd1 + inmon13 + m1];
                    xyz_sitei[1] = xyz_[fi_crd1 + inmon13 + nmon1 + m1];
                    xyz_sitei[2] = xyz_[fi_crd1 + inmon13 + 2 * nmon1 + m1];

                    for (size_t j = 0; j < ns2; j++) {
                        size_t jnmon2 = j * nmon2;
                        size_t jnmon23 = jnmon2 * 3;
                        double c6j = c6_long_range_[fi_sites2 + j * nmon2];
                        double c6, d6;
                        GetC6(mon_id_[fi_mon1], mon_id_[fi_mon2], i, j, c6, d6);
                        energy_pool[rank] +=
                            disp6(c6, d6, c6i, c6j, xyz_sitei.data(), xyz_.data() + fi_crd2, g1.data(),
                                  grad2_pool[rank].data(), phi_i, phi2_pool[rank].data(), nmon1, nmon2, m2init, nmon2,
                                  i, j, 1.0, do_grads_, cutoff_, ewald_alpha_, box_, box_inverse_, &virial_pool[rank]);
                    }
                    grad1_pool[rank][inmon13 + m1] += g1[0];
                    grad1_pool[rank][inmon13 + nmon1 + m1] += g1[1];
                    grad1_pool[rank][inmon13 + nmon12 + m1] += g1[2];
                    phi1_pool[rank][inmon1 + m1] += phi_i;
                }
            }

            // Compress data in Efq and phi
            for (size_t rank = 0; rank < nthreads; rank++) {
                size_t kend1 = grad1_pool[rank].size();
                size_t kend2 = grad2_pool[rank].size();
                for (size_t k = 0; k < kend1; k++) {
                    grad_[fi_crd1 + k] += grad1_pool[rank][k];
                }
                for (size_t k = 0; k < kend2; k++) {
                    grad_[fi_crd2 + k] += grad2_pool[rank][k];
                }
                kend1 = phi1_pool[rank].size();
                kend2 = phi2_pool[rank].size();
                for (size_t k = 0; k < kend1; k++) {
                    phi_[fi_sites1 + k] += phi1_pool[rank][k];
                }
                for (size_t k = 0; k < kend2; k++) {
                    phi_[fi_sites2 + k] += phi2_pool[rank][k];
                }
                for (size_t k = 0; k < 9; k++) {
                    virial_[k] += virial_pool[rank][k];
                }
                disp_energy_ += energy_pool[rank]; 
            }


            // Update first indexes
            fi_mon2 += nmon2;
            fi_sites2 += nmon2 * ns2;
            fi_crd2 += nmon2 * ns2 * 3;
        }
        // Update first indexes
        fi_mon1 += nmon1;
        fi_sites1 += nmon1 * ns1;
        fi_crd1 += nmon1 * ns1 * 3;
    }

    if (ewald_alpha_ > 0 && use_pbc_) {
        helpme::PMEInstance<double> pme_solver_;
        // Compute the reciprocal space terms, using PME
        double A = box_[0], B = box_[4], C = box_[8];
        int grid_A = pme_grid_density_ * A;
        int grid_B = pme_grid_density_ * B;
        int grid_C = pme_grid_density_ * C;
        pme_solver_.setup(6, ewald_alpha_, pme_spline_order_, grid_A, grid_B, grid_C, -1, 0);
        pme_solver_.setLatticeVectors(A, B, C, 90, 90, 90, PMEInstanceD::LatticeType::XAligned);
        // N.B. these do not make copies; they just wrap the memory with some metadata
        auto coords = helpme::Matrix<double>(sys_xyz_.data(), natoms_, 3);
        auto params = helpme::Matrix<double>(sys_c6_long_range_.data(), natoms_, 1);
        auto forces = helpme::Matrix<double>(sys_grad_.data(), natoms_, 3);
        std::vector<double> dummy_6vec(6,0.0);
        auto rec_virial = helpme::Matrix<double>(dummy_6vec.data(), 6, 1);
        std::fill(sys_grad_.begin(), sys_grad_.end(), 0);
        double rec_energy = pme_solver_.computeEFVRec(0, params, coords, forces, rec_virial);

        // get virial
        if (calc_virial_) {
            
            virial_[0] += *rec_virial[0];
            virial_[1] += *rec_virial[1];
            virial_[2] += *rec_virial[3];
            virial_[4] += *rec_virial[2];
            virial_[5] += *rec_virial[4];
            virial_[8] += *rec_virial[5];

            virial_[3] = virial_[1];
            virial_[6] = virial_[2];
            virial_[7] = virial_[5];

        }	

        // Resort forces from system order
        fi_mon = 0;
        fi_sites = 0;
        for (size_t mt = 0; mt < mon_type_count_.size(); mt++) {
            size_t ns = num_atoms_[fi_mon];
            size_t nmon = mon_type_count_[mt].second;
            for (size_t m = 0; m < nmon; m++) {
                size_t mns = m * ns;
                for (size_t i = 0; i < ns; i++) {
                    size_t inmon = i * nmon;
                    const double *result_ptr = forces[fi_sites + mns + i];
                    grad_[3 * fi_sites + 3 * inmon + 0 * nmon + m] -= result_ptr[0];
                    grad_[3 * fi_sites + 3 * inmon + 1 * nmon + m] -= result_ptr[1];
                    grad_[3 * fi_sites + 3 * inmon + 2 * nmon + m] -= result_ptr[2];
                }
            }
            fi_mon += nmon;
            fi_sites += nmon * ns;
        }
        // The Ewald self energy
        double prefac = std::pow(ewald_alpha_, 6) / 12.0;
        double self_energy = 0;
        for (const auto &c6 : sys_c6_long_range_) {
            self_energy += c6 * c6 * prefac;
        }
        disp_energy_ += rec_energy + self_energy;
    }
}

}  // namespace disp
