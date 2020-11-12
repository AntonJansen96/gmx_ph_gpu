/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
/*
 * constant_ph.h
 *
 */

#ifndef CONSTANT_PH_DATA_H
#define CONSTANT_PH_DATA_H

#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/random/threefry.h"
#include <stdio.h>
#include <config.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

/* Structure for one lambda group */
struct t_lambdarec
{
    /* x is for lambda coordinate */
    real                 x0;         /* initial lambda */
    real                 x;          /* lambda coordinate */
    real                 x_old;      /* lambda of previous step */
    real                 x_old_old;  /* lambda of previous step */
    real                 v;          /* velocity of lambda */
    real                 v_old;      /* velocity of previous step */
    real                 v_old_old;  /* velocity of previous step */
    real                 ekin;       /* kinetic energy of lambda */
    real                 T;          /* temperature */
    int                  tau;        /* tcouple time constant */
    real                 bar;        /* double well barrier */
    real                 dvdl;       /* dV/dlambda, force acting on lambda */
    real                 dvdl_old;   /* dV/dlambda, force acting on lambda */
    std::array<real, 15> lambda_dwp; /* double well potential array*/
    real                 dvdl_pot;   /* dvdl from environment */
    real                 dvdl_ph;    /* dvdl from ph */
    real                 dvdl_dwp;   /* dvdl from double well potential */
    real                 dvdl_ref;   /* dvdl from reference free energy */
};

/* Structure for lambda residues */
struct LambdaResidues
{
    std::string       residue_name;                   /* name of residue */
    int               residue_index  = 0;             /* reference to residues in cphmd_general*/
    int               residue_number = 0;             /* residue number */
    int               n_atoms        = 0;             /* number of atoms in lambda group */
    t_lambdarec       lambda;                         /* lambda structure */
    std::vector<real> chargeA;                        /* charge A */
    std::vector<real> chargeB;                        /* charge B */
    std::vector<int>  atoms;                          /* indexes of atoms of lambda group */
    FILE*             out = nullptr;                  /* output file for lambda dynamics data */
    int               chargeRestraintGroupNumber = 0; /* reference for charge restraint group */
    bool              isBufferResidue = false; /* if this residue is part of the buffer group */
};

struct dvdlCoeffecients
{
    std::vector<real> values_;
};

/* Structure for all lambda groups */
struct LambdaResidueCollection
{
    std::vector<LambdaResidues> residues_;
};

/* Structure for residue specific data */
struct cphmd_general
{
    real             ph_value              = 0.0; /* value of pH */
    int              nst_lambda            = 0;   /* output frequency*/
    real             tau_lambda            = 0.0; /* T couple tau for lambda */
    real             reference_temperature = 0.0;
    real             lambda_mass           = 0.0;
    std::vector<int> multiGroupLambdaStates; /* array for how many states in each multigroup lambda */
    std::vector<int> constrainedLambdaIndices;   /* indexes for lambdas one wants to constrain */
    int              n_buf           = 0;        /* number of (collective) buffers */
    real             m_buf           = 0.0;      /* mass of buffer*/
    real             lambda_tot_init = 0.0;      /* sum of lambdas for charge constraint */
    std::vector<std::string>      residue_names; /* array for residue names */
    std::vector<dvdlCoeffecients> dvdlCoefs;     /* array of dvdl values for each residue */
    std::vector<dvdlCoeffecients> edgeDvdlCoefs; /* array for dvdl_coeffs for multiple state edges */
    std::vector<real>             pKa_values;    /* array of pKa values */
    bool haveChargeConstraints     = false;      /* If we are using charge constraints */
    bool haveMultiStateConstraints = false;      /* If there are any multi state constraints */
    bool isCalibrationRun =
            false; /* If we are using this simulation to calibrate a set of dvdl parameters */
};
#endif
