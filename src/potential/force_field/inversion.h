#ifndef INVERSION_H
#define INVERSION_H

#include <math.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "topology.h"
#include "tools/custom_exceptions.h"

/**
 * @file inversion.h
 * @brief Function definitions for inversion class
 */

/**
 * @brief Inversion class constructs inversion objects and calculates the
 * potential energy or gradients
 *
 * The inversion class considers one functional form: harmonic
 * Its potentials is:
 * Harmonic: 0.5 * k * ( phi ijkn - phi 0 ) * ( phi ijkn - phi 0 )
 *
 * Further information on the potentials can be found at:
 * ftp://ftp.dl.ac.uk/ccp5/DL_POLY/DL_POLY_CLASSIC/DOCUMENTS/USRMAN2.19.pdf
 *
 * Parameter dictionary:
 * Harmonic:
 * linear_param[0] -> constant k
 * nonlinear_param[0] -> constant phi0
 *
 * The above parameter dictionary also shows how you should enter the parameters
 * in the connectivity file. Below is an example of what each number means
 * inversion     1         2        3       4      harm        1.0   2.0
 * ^             ^         ^        ^       ^      ^           ^      ^
 * topology   1st idx  2nd idx   3rd idx  4th idx  func form   k      phi0
 *
 */

class Inversion : public Topology {
   public:
    /**
     * Default constructor
     */
    Inversion();

    /**
     * @brief Overloaded constructor, setting the respective inversion_type
     * indexes, and functional_form.
     *
     * This constructor sets the field variables. It also puts the respective
     * number of 0's onto the vectors nonlinear_parameters_ and
     * linear_parameters_.
     * @param[in] topology This is the topology, which will be inversion
     * @param[in] indexes The indexes of the the atoms in the inversion angle
     * @param[in] functional_form The functional form that is used to evaluate
     *                            the energy
     */
    Inversion(std::string topology, std::vector<size_t> indexes, std::string functional_form);

    /**
     * Default Destructor.
     */
    ~Inversion();

    /**
     * @brief Calculates the potential energy using the nonlinear and linear
     *        parameters
     * @param[in] x Represents the inversion angle phi ijkn,
     * @return The potential energy for a given inversion angle phi ijkn
     */
    double GetEnergy(double x);

    /**
     * @brief gets the gradient for this particular topology for a given
     *        functional form
     * @param  x inversion angle phi ijkn
     * @return   the gradient
     *
     * Note. this is not the entire gradient. You still need to back propagate
     * this gradient to each of the individual x,y,z components
     */
    double GetTopologyGradient(double x);

    /**
     * @brief Checks if two inversions are the same by inspecting
     * all of the field variables.
     * @param[in] inversion Is the other inversion object that we are comparing
     * against
     * @return True if the inversion objects are the same. False otherwise
     */
    bool operator==(Inversion const &inversion) const;

   private:
};

#endif