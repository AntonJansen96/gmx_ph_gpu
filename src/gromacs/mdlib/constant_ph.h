/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef GMX_MDLIB_CONSTANT_PH_H
#define GMX_MDLIB_CONSTANT_PH_H

#include <memory>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"
#include "gromacs/random/threefry.h"

struct cphmd_general;
struct t_inputrec;
struct LambdaResidueCollection;
struct t_lambdarec;
struct t_mdatoms;
struct LambdaResidues;

/*! \libinternal
 * \brief This holds the constant pH data and provides methods for mdrun to interact with
 */
class ConstantPH
{
public:
    ConstantPH(const t_inputrec& ir, int natoms, const std::string& inputFileName);

    ~ConstantPH();

    //! Return the buffer to add electrostatic potential contributions to, to beindexed by local atom index
    gmx::ArrayRef<real> potential() { return potential_; }

    //! Returns the lambda atom indices
    gmx::ArrayRef<const int> lambdaAtoms() const { return lambdaAtoms_; }

    //! Sets mdatoms->chargeA for particles coupled to lambda's
    void setLambdaCharges(t_mdatoms* mdatoms);

    //! Update the lambda variables using the computed potential
    //
    // Returns the change in kinetic energy due to T-coupling
    //
    // TODO: Remove argument ir
    real updateLambdas(const t_inputrec& ir, double t, int64_t step);

    //! Need ugly access like this to populate data from checkpoint.
    LambdaResidueCollection* lambdaResidueHandle() { return ml_.get(); }

private:
    std::unique_ptr<cphmd_general> cphmd_gen_;
    // TODO: Replace these plain pointers
    std::unique_ptr<LambdaResidueCollection> ml_;
    std::vector<int>                         lambdaAtoms_;
    std::vector<int>                         lambdaAtomsIndex_;
    std::vector<real>                        chargeA_;
    std::vector<real>                        chargeB_;
    std::vector<real>  dvdl_; // dVdl sum for each lambda group (size = number of lambdas)
    std::vector<real>  potential_;
    const std::string& inputFileName_;
    int                eLambdaThermostat_ = 0; // is equal to eLambdaTcVRESCALE, ugly but effective.
    bool               useChargeConstraints_    = false;
    bool               useMultiStateConstraits_ = false;
    bool               isCalibrationRun_        = false;
};

constexpr int MAX_N_DVDL_COEFFS = 10;
void          update_lambda(const t_inputrec&         ir,
                            int64_t                   step,
                            const cphmd_general&      cphmd_gen,
                            LambdaResidues*           ml,
                            gmx::ArrayRef<const real> pot);


#endif
