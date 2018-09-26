/*******************************************************************************
This file contains a class that allows to compute the permanent electric field,
the dipole field, and the gradients for a given pair of sites over all the 
monomers. In the same calculation it modifies fields, electric fields and 
gradients of both sites, returning by reference the values of site i from 
mon 1 and modifying in situ the corresponding ones for site 2. This enables 
vectorization.

!!! IMPORTANT. READ. REALLY, READ. !!!
This functions assume that the data in all the pointers is already CONTIGIOUS
The order assumed is the following:

- XYZ, MU, EFQ, EFD, GRD
FIXME fix this description. This is false. Look at the reorganization in electrostatics.cpp
[BLOCK1][BLOCK2]...[BLOCKN]
where each block is one monomer type organized in the way:
x.1.1 x.1.2 ... x.2.1. x.2.2 ... x.A.B ... y.1.1 y.1.2 ... z.1.1 ... z.A.B
where x.1.2 will be read as coordinate x of site 1 of monomer 2, for B 
monomers with A sites 

- PHI, CHG
[BLOCK1][BLOCK2]...[BLOCKN]
where each block is one monomer type organized in the way:
phi.1.1 phi.1.2 ... phi.2.1 ... phi.A.B 
where phi.1.2 will be read as potential on site 1 of monomer 2, for B 
monomers with A sites

!!! END IMPORTANT !!!

The equations are taken from the manuscript by Masia 
dx.doi.org/10.1063/1.3511713 , J. Chem. Phys. 133, 234101 (2010)

NOTE: The screening functions in this paper are for TTM3, while here we are
using TTM4. The exponents change from 3 to 4 but everything can be trivially 
derived.
*******************************************************************************/
#ifndef ELECTROSTATIC_TENSORS_H
#define ELECTROSTATIC_TENSORS_H

#include <vector>
#include <cmath>
#include <cstddef>

#include "potential/electrostatics/gammq.h" 

namespace elec {

////////////////////////////////////////////////////////////////////////////////

class ElectroTensor {
 public:
  ElectroTensor() {maxnmon = 0;}
  ElectroTensor(size_t n);
  ~ElectroTensor() {}

////////////////////////////////////////////////////////////////////////////////
// Get tensor alpha ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  // Computes the electric field for a pair of sites for a number of monomers
  // # = mon2_index_end - mon2_index_start when A=polfac[i] * polfac[j] > 0

  void CalcT0AndT1(
    double * xyz1, double * xyz2,     // Coordinates of mon type 1 and 2
    size_t mon1_index,                // Mon 1 index
    size_t mon2_index_start,          // Mon 2 initial index
    size_t mon2_index_end,            // Mon 2 final index
    size_t nmon1, size_t nmon2,       // # monomers of types 1 and 2
    size_t site_i, size_t site_j,     // Site # i of mon1 and # j of mon 2
    double Ai,                        // (polfac[i] * polfac[j]) inverted
    double Asqsqi,                     // (polfac[i] * polfac[j])^4
    double aCC, double aCC1_4,        // Thole damping aCC and aCC^(0.25)
    double g34,                       // Gamma ln function. Is a constant.
    // The following pointers correspond to the same array, but they point to
    // the part of the array that contains information of mons 1 and 2
    double * ts0_mon1, // Size n*n    // Electrostatic tensor Tij 
    double * ts0_mon2, // Size n*n    // Electrostatic tensor Tij
    double * ts1_mon1, // Size 3n*n   // Electrostatic tensor Tij_alpha for Mon 1
    double * ts1_mon2);// Size 3n*n   // Electrostatic tensor Tij_alpha for Mon 2

////////////////////////////////////////////////////////////////////////////////
// Get Tensor alpha beta ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  // Computes the dipole field for a pair of sites for a number of monomers
  // # = mon2_index_end - mon2_index_start when A=polfac[i] * polfac[j] > 0

  void CalcT1AndT2(
    double * xyz1, double * xyz2,     // Coordinates of mon type 1 and 2
    size_t mon1_index,                // Mon 1 index
    size_t mon2_index_start,          // Mon 2 initial index
    size_t mon2_index_end,            // Mon 2 final index
    size_t nmon1, size_t nmon2,       // # monomers of types 1 and 2
    size_t site_i, size_t site_j,     // Site # i of mon1 and # j of mon 2
    double Asqsqi,                     // (polfac[i] * polfac[j])^4
    double aDD,                       // Thole damping aDD (dipole - dipole)
    size_t nsites,
    // The following pointers correspond to the same array, but they point to
    // the part of the array that contains information of mons 1 and 2
    double * ts1_mon1,  // Size 3n*n  // Electrostatic tensor Tij_alpha for Mon 1
    double * ts1_mon2,  // Size 3n*n  // Electrostatic tensor Tij_alpha for Mon 2
    double * ts2); // Size 3n*3n // Electrostatic tensor Tij_alpha_beta

////////////////////////////////////////////////////////////////////////////////

 private:

  // Maximum number of monomers that can be in system 
  // For now, is set to the largest number. Later can be set to
  // the maximum number of monomers we will evaluate at once.
  size_t maxnmon; 
};

} // namespace elec
#endif
