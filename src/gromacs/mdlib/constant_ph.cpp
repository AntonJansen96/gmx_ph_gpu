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
/*! \file
 * \brief
 * Implements constant pH lambda dynamics code.
 *
 * \author Noora Aho <noora.s.aho@jyu.fi>
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 *
 */

#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdlib/constant_ph.h"
#include "gromacs/mdlib/constant_ph_data.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/random/seed.h"
#include "gromacs/random/normaldistribution.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace
{

constexpr real s             = 0.3;
constexpr real one           = 1.0;
constexpr real two           = 2.0;
constexpr real ten           = 10.0;
constexpr real half          = 0.5;
constexpr real oneTenth      = 0.1;
constexpr real oneFifth      = 0.2;
constexpr real oneTwentieth  = 0.05;
constexpr real oneHundredth  = 0.01;
constexpr real oneThousendth = 0.001;
constexpr real zero          = 0.0;

/* Initialize double well barrier potential
   Previous work taken from constant pH version 5.1
 */
void init_lambda_dwp(gmx::ArrayRef<real> lambda_dwp, const real barrier)
{

    /* fields of array lambda_dwp are filled 08.08.2014 */
    constexpr real eps      = 0.005;
    int            iter     = 0;
    constexpr int  max_iter = 10000;

    const real sig0  = oneFifth;
    real       sig   = sig0;
    real       x0    = zero;
    const real vmin0 = ten;

    /* Note that the height of outer walls will be f+barrier/2 */
    // For f=1000 follows
    constexpr real f   = 1000.0;
    constexpr real e1  = 2.185;
    constexpr real e10 = 1.645;

    lambda_dwp[0] = s;
    lambda_dwp[1] = f;
    lambda_dwp[2] = e1;
    lambda_dwp[3] = e10;

    const real st = (e1 - e10) / (two * sig0);
    const real m  = two * sig0 * (two * e1 - e10) / (e1 - e10);

    lambda_dwp[4] = st;
    lambda_dwp[5] = m;

    const real flagfor = barrier - one;

    for (real iii = flagfor; iii <= flagfor; iii = iii + 1)
    {
        const real dg   = one + iii;
        const real sig0 = 0.02;
        sig             = sig0;
        real b          = -oneTenth;
        real a          = oneTwentieth;
        real k          = dg / two;
        real d          = dg / two;

        /*  v=-k*(expf(-((x-1.0-b)*(x-1.0-b)/(2*a*a)))+expf(-(x+b)*(x+b)/(2*a*a)))+
         *       d*expf(-(x-0.5)*(x-0.5)/(2*s*s))+
         *            f*0.5*((1-std::erf(st*(x+m)))+(1.0+std::erf(st*(x-1.0-m))));
         *                this is the barrier function */

        /* correct for the fact that the two minima are shallower than k=dg/2 due to the tails of
         * the central gaussian and the erf functions (basically bring the minima to the same value
         * of the two side gaussians before the central gaussian was added) */

        real vmin = vmin0;
        for (real x = -oneTenth; x <= oneFifth; x = x + oneThousendth)
        {
            const real v =
                    -k
                            * (std::exp(-((x - one - b) * (x - one - b) / (two * a * a)))
                               + std::exp(-(x + b) * (x + b) / (two * a * a)))
                    + d * std::exp(-(x - half) * (x - half) / (two * s * s))
                    + f * -half * ((one - std::erf(st * (x + m))) + (one + std::erf(st * (x - one - m))));
            if (v < vmin)
            {
                vmin = v;
            }
        }
        k = k + dg / two + vmin;

        /* adjust location minima and width boltzmann distribution therein to the target values (0,1 and sig0) */

        bool ready = false;
        while (!ready) /*while which correspond to repeat begin*/

        {
            real vmin = vmin0;
            for (real x = -oneTenth; x < oneFifth; x = x + oneThousendth)
            {
                const real v =
                        -k
                                * (std::exp(-((x - 1.0 - b) * (x - 1.0 - b) / (2 * a * a)))
                                   + std::exp(-(x + b) * (x + b) / (2 * a * a)))
                        + d * std::exp(-(x - 0.5) * (x - 0.5) / (2 * s * s))
                        + f * 0.5 * ((1 - std::erf(st * (x + m))) + (1.0 + std::erf(st * (x - 1.0 - m))));
                if (v < vmin)
                {
                    vmin = v;
                }
            }
            real totp = 0;
            {
                k         = k + dg / two + vmin;
                real x    = half;
                real totx = 0;
                while (x > -oneFifth)
                {
                    const real v = -k
                                           * (std::exp(-((x - 1.0 - b) * (x - 1.0 - b) / (2 * a * a)))
                                              + std::exp(-(x + b) * (x + b) / (2 * a * a)))
                                   + d * std::exp(-(x - 0.5) * (x - 0.5) / (2 * s * s))
                                   + f * 0.5
                                             * ((1 - std::erf(st * (x + m)))
                                                + (1.0 + std::erf(st * (x - 1.0 - m))));
                    if (v <= 0)
                    {
                        totp = totp + std::exp(-v);
                        totx = totx + x * std::exp(-v);
                    }
                    x = x - oneThousendth;
                }
                x0 = totx / totp;
            }

            {
                real x     = half;
                real totxp = 0;
                while (x > -oneFifth)
                {
                    const real v = -k
                                           * (std::exp(-((x - 1.0 - b) * (x - 1.0 - b) / (2 * a * a)))
                                              + std::exp(-(x + b) * (x + b) / (2 * a * a)))
                                   + d * std::exp(-(x - 0.5) * (x - 0.5) / (2 * s * s))
                                   + f * 0.5
                                             * ((1 - std::erf(st * (x + m)))
                                                + (1.0 + std::erf(st * (x - 1.0 - m))));
                    if (v <= 0)
                    {
                        totxp = totxp + ((x - x0) * (x - x0) * std::exp(-v));
                    }
                    x = x - oneThousendth;
                }
                sig = std::sqrt(totxp / totp);
                b   = b + oneHundredth * x0;
                a   = a / (1.0 + oneHundredth * (sig - sig0) / sig0);
            }

            if ((std::abs(x0) <= eps) && (std::abs(sig - sig0) / sig0 <= eps))
            {
                ready         = true;
                lambda_dwp[6] = a;
                lambda_dwp[7] = b;
                lambda_dwp[8] = k;
                lambda_dwp[9] = d;
            }

            if (iii < -(half + oneTwentieth))
            {
                k             = 0;
                d             = 0;
                ready         = true;
                lambda_dwp[6] = a;
                lambda_dwp[7] = b;
                lambda_dwp[8] = k;
                lambda_dwp[9] = d;
                if (barrier != 0)
                {
                    fprintf(stderr, "Warning: Barrier changed from %f to zero\n", barrier);
                }
                fprintf(stderr, "Parameters double well potential:\n");
                fprintf(stderr, "a = %f b = %f k = %f d = %f dg = %f \n", a, b, k, d, dg);
                fprintf(stderr, "m = %f st = %f \n", m, st);
            }

            iter++;

            if (iter > max_iter)
            {
                ready         = true;
                lambda_dwp[6] = a;
                lambda_dwp[7] = b;
                lambda_dwp[8] = k;
                lambda_dwp[9] = d;
                fprintf(stderr,
                        "Warning: double well potential did not converge to tolerance limit eps=%f "
                        "in max iter=%d\n",
                        eps, max_iter);
                fprintf(stderr, "Parameters double well potential:\n");
                fprintf(stderr, "a = %f b = %f k = %f d = %f dg = %f \n", a, b, k, d, dg);
                fprintf(stderr, "m = %f st = %f \n", m, st);
            }

        } /* endwhile of repeat begin */

    } /* end for */
}

real gaussdist(gmx::DefaultRandomEngine* rng, real sigma)
{
    constexpr real                     two = 2.0;
    real                               r   = two;
    gmx::UniformRealDistribution<real> uniformDist;
    real                               x = 0;
    do
    {
        x            = two * uniformDist(*rng) - 1.0;
        const real y = 2.0 * uniformDist(*rng) - 1.0;
        r            = x * x + y * y;
    } while (r > 1.0 || r == 0.0);
    r = x * std::sqrt(-two * std::log(r) / r);
    r = r * sigma;
    return (r);
}

using namespace std;

/*  TODO: Maybe add more input checking and errors for input */
void init_constantph(LambdaResidueCollection* ml,
                     cphmd_general*           cphmd_gen,
                     const std::string&       inputFileName,
                     bool                     useChargeConstraints,
                     bool                     useMultiStateConstraits,
                     bool                     isCalibrationRun)
{
    fprintf(stderr, "\nConstant pH MD initialization: \n\n");

    int                   numLambdaGroups   = 0;
    int                   numLambdaResidues = 0;
    ifstream              myfile(inputFileName);
    std::array<string, 3> temp;
    string                line;
    stringstream          ss;

    getline(myfile, line);
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->ph_value;
    cout << temp[0] << " " << cphmd_gen->ph_value << "\n";
    ss.clear();

    getline(myfile, line);
    ss << line;
    ss >> temp[0] >> temp[1] >> numLambdaResidues;
    cout << temp[0] << " " << numLambdaResidues << "\n";
    ss.clear();

    getline(myfile, line);
    ss << line;
    ss >> temp[0] >> temp[1] >> numLambdaGroups;
    cout << temp[0] << " " << numLambdaGroups << "\n";
    ss.clear();

    getline(myfile, line);

    getline(myfile, line);
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->lambda_mass;
    cout << temp[0] << " " << cphmd_gen->lambda_mass << "\n";
    ss.clear();

    getline(myfile, line);
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->reference_temperature;
    cout << temp[0] << " " << cphmd_gen->reference_temperature << "\n";
    ss.clear();

    getline(myfile, line);
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->tau_lambda;
    cout << temp[0] << " " << cphmd_gen->tau_lambda << "\n";
    ss.clear();

    getline(myfile, line);
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->nst_lambda;
    cout << temp[0] << " " << cphmd_gen->nst_lambda << "\n";
    ss.clear();

    cphmd_gen->haveMultiStateConstraints = useMultiStateConstraits;

    int sumOfLambdaStates = 0;
    if (useMultiStateConstraits)
    {
        fprintf(stderr, "\nUsing multi state constraints\n\n");

        int numMultiStateConstraintGroups = 0;
        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1] >> numMultiStateConstraintGroups;
        cout << temp[0] << " " << numMultiStateConstraintGroups << "\n\n";
        ss.clear();

        cphmd_gen->multiGroupLambdaStates.resize(numMultiStateConstraintGroups);

        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1];
        cout << temp[0] << " " << temp[1] << " ";
        for (auto& multiGroupLambdaState : cphmd_gen->multiGroupLambdaStates)
        {
            ss >> multiGroupLambdaState;
            cout << multiGroupLambdaState << " ";
            sumOfLambdaStates += multiGroupLambdaState;
        }

        ss.clear();
        cout << "\n";
    }

    cphmd_gen->haveChargeConstraints = useChargeConstraints;
    if (useChargeConstraints)
    {
        fprintf(stderr, "\nUsing charge constraints\n\n");
        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1] >> cphmd_gen->n_buf;
        cout << temp[0] << " " << cphmd_gen->n_buf << "\n";
        ss.clear();

        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1] >> cphmd_gen->m_buf;
        cout << temp[0] << " " << cphmd_gen->m_buf << "\n";
        ss.clear();

        int numConstrainedLambdas =
                numLambdaResidues + cphmd_gen->multiGroupLambdaStates.size() - sumOfLambdaStates;
        cphmd_gen->constrainedLambdaIndices.resize(numConstrainedLambdas);

        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1];
        cout << temp[0] << " " << temp[1] << " ";
        for (auto& constrainedLambdaIndex : cphmd_gen->constrainedLambdaIndices)
        {
            ss >> constrainedLambdaIndex;
            cout << constrainedLambdaIndex << " ";
        }
        ss.clear();

        cout << "\n\n";
    }


    cphmd_gen->residue_names.resize(numLambdaResidues);
    cphmd_gen->pKa_values.resize(numLambdaResidues);
    cphmd_gen->dvdlCoefs.resize(numLambdaResidues);

    /* Loop over residue inputs */
    for (int j = 0; j < numLambdaResidues; j++)
    {
        getline(myfile, line); // empty line

        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1] >> cphmd_gen->residue_names[j];
        cout << temp[0] << " " << cphmd_gen->residue_names[j] << "\n";
        ss.clear();

        int numCoeff = 0;
        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1] >> numCoeff;
        cout << temp[0] << " " << numCoeff << "\n";
        ss.clear();
        cphmd_gen->dvdlCoefs[j].values_.resize(numCoeff);

        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1];
        cout << temp[0] << " " << temp[1] << " ";
        for (auto& coeff : cphmd_gen->dvdlCoefs[j].values_)
        {
            ss >> coeff;
            cout << coeff << " ";
        }
        ss.clear();

        int numCoeffEdge = 0;
        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1] >> numCoeffEdge;
        cout << temp[0] << " " << numCoeffEdge << "\n";
        ss.clear();
        cphmd_gen->edgeDvdlCoefs[j].values_.resize(numCoeffEdge);

        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1];
        cout << temp[0] << " " << temp[1] << " ";
        for (auto& coeff : cphmd_gen->edgeDvdlCoefs[j].values_)
        {
            ss >> coeff;
            cout << coeff << " ";
        }
        ss.clear();
        cout << "\n";

        cout << "\n";

        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1] >> cphmd_gen->pKa_values[j];
        cout << temp[0] << " " << cphmd_gen->pKa_values[j] << "\n\n";
        ss.clear();
    }

    /* Loop over lambda groups ml->next */
    int pos                    = 0;
    cphmd_gen->lambda_tot_init = 0;

    for (int j = 0; j < numLambdaGroups; j++)
    {
        LambdaResidues newResidue;

        getline(myfile, line); // empty line

        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1] >> newResidue.residue_name;
        cout << "Residue name set through " << temp[0] << " " << newResidue.residue_name << "\n";
        ss.clear();

        // check which residue in cphmd_gen the lambda group relates to
        for (int k = 0; k < numLambdaResidues; k++)
        {
            if (newResidue.residue_name == cphmd_gen->residue_names[k])
            {
                newResidue.residue_index = k;
                break;
            }
        }

        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1] >> newResidue.residue_number;
        cout << temp[0] << " " << newResidue.residue_number << "\n";
        ss.clear();

        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1] >> newResidue.chargeRestraintGroupNumber;
        cout << temp[0] << " " << newResidue.chargeRestraintGroupNumber << "\n";
        ss.clear();

        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1] >> newResidue.lambda.x0;
        cout << temp[0] << " " << newResidue.lambda.x0 << "\n";
        ss.clear();

        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1] >> newResidue.lambda.bar;
        cout << temp[0] << " " << newResidue.lambda.bar << "\n";
        ss.clear();

        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1] >> newResidue.n_atoms;
        cout << temp[0] << " " << newResidue.n_atoms << "\n\n";
        ss.clear();

        newResidue.atoms.resize(newResidue.n_atoms);

        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1];
        for (auto atom = newResidue.atoms.begin(); atom != newResidue.atoms.end(); atom++)
        {
            ss >> *atom;
            *atom -= 1; // to shift index starting from 0
        }
        ss.clear();
        // Information about buffer residue or not.
        getline(myfile, line);
        ss << line;
        ss >> temp[0] >> temp[1] >> newResidue.isBufferResidue;
        cout << temp[0] << " " << newResidue.isBufferResidue << "\n";
        ss.clear();
        // empty line
        getline(myfile, line);

        // get charges of atoms in lambda group
        newResidue.chargeA.resize(newResidue.n_atoms);
        newResidue.chargeB.resize(newResidue.n_atoms);
        for (int kk = 0; kk < gmx::ssize(newResidue.atoms); kk++)
        {
            getline(myfile, line);
            ss << line;
            ss >> temp[0] >> newResidue.chargeA[kk] >> newResidue.chargeB[kk];
            cout << temp[0] << " " << newResidue.chargeA[kk] << " " << newResidue.chargeB[kk] << "\n";
            ss.clear();
        }
        cout << "\n";

        pos++;
        // open output file for this lambda
        std::string lambdaFileName = gmx::formatString("lambda_%d.dat", (pos));
        cout << lambdaFileName << "\n\n";
        newResidue.out = fopen(lambdaFileName.c_str(), "w");
        fprintf(newResidue.out, "# Name %s, Residue Number %d, id %d\n",
                newResidue.residue_name.c_str(), newResidue.residue_number, pos);
        fprintf(newResidue.out,
                "# Time, Lambda coordinate value, Lambda dvdl, Lambda Temperature\n");

        // initialize lambda structure
        newResidue.lambda.x    = newResidue.lambda.x0;
        newResidue.lambda.T    = cphmd_gen->reference_temperature;
        newResidue.lambda.tau  = cphmd_gen->tau_lambda;
        newResidue.lambda.dvdl = 0.0;

        if (!isCalibrationRun)
        {
            // initial velocity of lambda particle
            gmx::UniformRealDistribution<real> uniformDist;
            static int                         start_seed = 0;
            start_seed                                    = static_cast<int>(gmx::makeRandomSeed());
            gmx::DefaultRandomEngine Rng(start_seed);
            double sigma = sqrt(newResidue.lambda.T * BOLTZ / cphmd_gen->lambda_mass);
            if (newResidue.isBufferResidue)
            {
                newResidue.lambda.v     = 0.0;
                newResidue.lambda.v_old = 0.0;
            }
            else
            {
                newResidue.lambda.v = gaussdist(&Rng, sigma); /* random start velocity to lambda */
                newResidue.lambda.v_old = 0.0;
            }
        }
        else
        {
            newResidue.lambda.v     = 0.0;
            newResidue.lambda.v_old = 0.0;
        }

        // initialize double well potential
        init_lambda_dwp(newResidue.lambda.lambda_dwp, newResidue.lambda.bar);
        // initial sum of lambdas for charge constraint + buffer mass
        for (const auto& constrainedLambda : cphmd_gen->constrainedLambdaIndices)
        {
            if (pos == constrainedLambda)
            {
                if (newResidue.isBufferResidue)
                {
                    cphmd_gen->lambda_tot_init =
                            cphmd_gen->lambda_tot_init
                            + static_cast<real>(cphmd_gen->n_buf) * newResidue.lambda.x0;
                }
                else
                {
                    cphmd_gen->lambda_tot_init = cphmd_gen->lambda_tot_init + newResidue.lambda.x0;
                }
            }
        }

        ml->residues_.push_back(newResidue);
    }

    if (isCalibrationRun)
    {
        fprintf(stderr,
                "Current run will calibrate a set of dvdl parameters and not update the lambda "
                "coordinate\n");
    }
    cphmd_gen->isCalibrationRun = isCalibrationRun;

    fprintf(stderr, "lambda_tot_init = %f\n", cphmd_gen->lambda_tot_init);
    fprintf(stderr, "\n");
}

void compute_forces(const cphmd_general&                cphmd_gen,
                    LambdaResidues*                     currentLambdaResidue,
                    gmx::ArrayRef<const LambdaResidues> lambdaResidues,
                    real                                dvdl,
                    const int                           currentLamba)
{
    // real  dvdl = 0.0;
    real corr_ph   = 0.0;
    real corr_dvdl = 0.0;
    real corr_dwp  = 0.0;

    real                      ref_pka = cphmd_gen.pKa_values[currentLambdaResidue->residue_index]; // reference pKa for current lambda group
    gmx::ArrayRef<const real> dvdlCoefs = cphmd_gen.dvdlCoefs[currentLambdaResidue->residue_index].values_;

    int numberOfConstraintGroups = gmx::ssize(cphmd_gen.multiGroupLambdaStates);
    // effect of bond breaking (free energy dvdl) -> multi state or simple two state
    const int lambdaChargeRestraintIndex = currentLambdaResidue->chargeRestraintGroupNumber - 1;
    if (cphmd_gen.haveMultiStateConstraints
        && currentLambdaResidue->chargeRestraintGroupNumber <= numberOfConstraintGroups
        && cphmd_gen.multiGroupLambdaStates[lambdaChargeRestraintIndex] > 2)
    {
        std::vector<real> lambdaArray;
        // save two other lambda coordinates for this lambda dV/dl
        for (int i = 0; i < gmx::ssize(lambdaResidues); i++)
        {
            if (i == currentLamba)
            {
                continue;
            }
            lambdaArray.emplace_back(lambdaResidues[i].lambda.x);
        }

        gmx::ArrayRef<const real> edgeDvdlCoefs =
                cphmd_gen.edgeDvdlCoefs[currentLambdaResidue->residue_index].values_;
        ;
        // check which dV/dl to pick, which side of constrained triangle
        if ((lambdaArray[currentLamba] + lambdaArray[0]) > (lambdaArray[currentLamba] + lambdaArray[1]))
        {
            int i = 0;
            for (const auto& coeff : dvdlCoefs)
            {
                corr_dvdl += -1.0 * coeff * pow(currentLambdaResidue->lambda.x, i);
                i++;
            }
        }
        else
        {
            int i = 0;
            for (const auto& coeff : edgeDvdlCoefs)
            {
                corr_dvdl += -1.0 * coeff * pow(currentLambdaResidue->lambda.x, i);
                i++;
            }
        }
    }

    // dvdl array now has the force so no need to loop over atoms (lambdaAtoms array can be ordered)

    // effect of external pH
    corr_ph += BOLTZ * cphmd_gen.reference_temperature * std::log(ten) * (ref_pka - cphmd_gen.ph_value);

    // effect of double barrier potential - taken from previous constant ph version 5.1
    real xx = currentLambdaResidue->lambda.x;
    real s  = currentLambdaResidue->lambda.lambda_dwp[0];
    real f  = currentLambdaResidue->lambda.lambda_dwp[1];
    real st = currentLambdaResidue->lambda.lambda_dwp[4];
    real m  = currentLambdaResidue->lambda.lambda_dwp[5];
    real a  = currentLambdaResidue->lambda.lambda_dwp[6];
    real b  = currentLambdaResidue->lambda.lambda_dwp[7];
    real k  = currentLambdaResidue->lambda.lambda_dwp[8];
    real d  = currentLambdaResidue->lambda.lambda_dwp[9];

    corr_dwp = -k
                       * ((-1.0 * (xx - 1.0 - b)
                           * expf(-half * (xx - 1.0 - b) * (xx - 1.0 - b) / (a * a)) / (a * a))
                          + (-1.0 * (xx + b) * expf(-half * (xx + b) * (xx + b) / (a * a)) / (a * a)))
               - d * (xx - half) * expf(-(xx - half) * (xx - half) / (2 * s * s)) / (s * s)
               + f * half
                         * (2 * st * expf(-(st * st) * (xx - 1.0 - m) * (xx - 1.0 - m)) / std::sqrt(M_PI)
                            - 2 * st * expf(-(st * st) * (xx + m) * (xx + m)) / std::sqrt(M_PI));

    // for dvdl output
    currentLambdaResidue->lambda.dvdl_dwp = corr_dwp;
    currentLambdaResidue->lambda.dvdl_ref = corr_dvdl;
    currentLambdaResidue->lambda.dvdl_ph  = corr_ph;
    currentLambdaResidue->lambda.dvdl_pot = dvdl;

    // If mass zero then don't update -> reference free energy calculation
    if (!cphmd_gen.isCalibrationRun)
    {
        dvdl += corr_ph + corr_dvdl + corr_dwp;
    }
    currentLambdaResidue->lambda.dvdl = dvdl;
}

/*
 *  Update lambda and velocity using stochastic (Langevin) dynamics
 *  - modified from gromacs langevin function doSDUpdateGeneral() in update.cpp
 *  - phases 1,2 take into account the charge constraint
 */
template<int phase>
void updateLambdaLD(const t_inputrec& ir, const cphmd_general& cphmd_gen, LambdaResidues* ml, const real gamma)
{

    const real T_ref = cphmd_gen.reference_temperature;
    if (!cphmd_gen.isCalibrationRun)
    {
        GMX_RELEASE_ASSERT(cphmd_gen.lambda_mass != 0.0,
                           "Mass can't be zero for non calibration runs");
        real inverseMass = 1.0 / cphmd_gen.lambda_mass;

        // compare to updateType == SDUpdate::ForcesOnly in update.cpp
        if (phase == 1)
        {
            ml->lambda.x_old_old = ml->lambda.x;
            ml->lambda.v_old_old = ml->lambda.v;

            ml->lambda.v = ml->lambda.v + inverseMass * (-1.0) * ml->lambda.dvdl * ir.delta_t;
            ml->lambda.x = ml->lambda.x + ml->lambda.v * ir.delta_t;

            ml->lambda.x_old = ml->lambda.x;
            ml->lambda.v_old = ml->lambda.v;
        }

        // compare to updateType == SDUpdate::FrictionAndNoiseOnly in update.cpp
        if (phase == 2)
        {
            int seed = static_cast<int>(gmx::makeRandomSeed());
            // copied from update.cpp
            gmx::ThreeFry2x64<0> rng(seed, gmx::RandomDomain::UpdateCoordinates);
            gmx::TabulatedNormalDistribution<real, 14> dist;

            const real kT = BOLTZ * T_ref;

            real inverseSqrMass = std::sqrt(inverseMass);

            const real f = 1.0 - std::exp(-1.0 * gamma * ir.delta_t);

            ml->lambda.v_old_old = ml->lambda.v;

            const real deltav = (1.0 - f) * ml->lambda.v
                                + inverseSqrMass * std::sqrt(kT * f * (2.0 - f)) * dist(rng);
            ml->lambda.x = ml->lambda.x + half * (deltav - ml->lambda.v) * ir.delta_t;

            ml->lambda.x_old_old = ml->lambda.x;

            // ml->lambda->v_old = 0.5*(deltav - ml->lambda->v);
            ml->lambda.v_old = deltav;
            ml->lambda.x_old = ml->lambda.x;
        }

        // in update.cpp the "normal" update without constraint, add friction and noise
        // here called if charge constraint is not on
        if (phase == 3)
        {
            int seed = static_cast<int>(gmx::makeRandomSeed());
            // copied from update.cpp
            gmx::ThreeFry2x64<0> rng(seed, gmx::RandomDomain::UpdateCoordinates);
            gmx::TabulatedNormalDistribution<real, 14> dist;

            const real kT = BOLTZ * T_ref;

            const real inverseSqrMass = std::sqrt(inverseMass);

            const real f = 1.0 - std::exp(-1.0 * gamma * ir.delta_t);

            ml->lambda.v      = ml->lambda.v + inverseMass * (-1.0) * ml->lambda.dvdl * ir.delta_t;
            const real deltav = -1.0 * f * ml->lambda.v
                                + inverseSqrMass * std::sqrt(kT * f * (2.0 - f)) * dist(rng);

            ml->lambda.x = ml->lambda.x + (ml->lambda.v + half * deltav) * ir.delta_t;
            ml->lambda.v = ml->lambda.v + deltav;

            ml->lambda.ekin = half * cphmd_gen.lambda_mass * ml->lambda.v * ml->lambda.v;
            ml->lambda.T    = two * (ml->lambda.ekin) / BOLTZ;
        }
    }
}

/*
 * Apply constraints for multiple states compound and total charge of the system
 * - now we assume that v-rescale thermostat is used (langevin support maybe later)
 *
 * Charge constraint with collective buffer from:
 *
 * Charge-Neutral Constant pH Molecular Dynamics Simulations Using a Parsimonious Proton Buffer
 * Serena Donnini, R.Thomas Ullmann, Gerrit Groenhof and Helmut Grubmüller
 * J. Chem. Theory Comput. 2016, 12, 1040−1051
 *
 */
void do_constraints(const cphmd_general& cphmd_gen, gmx::ArrayRef<LambdaResidues> lambdaResidues, real dt)
{
    // (1) multiple states current lambda

    constexpr real    multiLambdaInitial       = 1.0;
    const int         numberOfMultiStateGroups = cphmd_gen.multiGroupLambdaStates.size();
    const int         numConstrainedLambdas    = cphmd_gen.constrainedLambdaIndices.size();
    std::vector<real> deviation(numberOfMultiStateGroups, 0);
    real              max_deviation = 0.0;

    if (cphmd_gen.haveMultiStateConstraints)
    {
        std::vector<real> multiLambdaGroupState(numberOfMultiStateGroups, 0);
        // compute current sum of lambdas for all multi state groups
        for (int group = 1; group <= numberOfMultiStateGroups; group++)
        {
            for (const auto& lambdaResidue : lambdaResidues)
            {
                // check if lambda groups belongs to certain multiple states group
                if (lambdaResidue.chargeRestraintGroupNumber == group)
                {
                    multiLambdaGroupState[group - 1] += lambdaResidue.lambda.x;
                }
            }
        }

        for (int i = 0; i < numberOfMultiStateGroups; i++)
        {
            deviation[i] = abs(multiLambdaGroupState[i] - multiLambdaInitial);
        }
        max_deviation = *std::max_element(deviation.begin(), deviation.begin() + numberOfMultiStateGroups
                                                                     - 1); // largest deviation
    }

    // (2) charge constraint current total lambda

    const real lambdaChargeConstraintInitial = cphmd_gen.lambda_tot_init;
    real       lambdaChargeConstraint        = lambdaChargeConstraintInitial;

    if (cphmd_gen.haveChargeConstraints)
    {
        lambdaChargeConstraint = 0;
        // compute current sum of all lambdas and collective buffer
        for (int i = 0; i < gmx::ssize(lambdaResidues); i++)
        {
            const auto& lambdaResidue = lambdaResidues[i];
            for (int j = 0; j < numConstrainedLambdas; j++)
            {
                if (i == cphmd_gen.constrainedLambdaIndices[j] - 1)
                {
                    if (lambdaResidue.isBufferResidue)
                    {
                        lambdaChargeConstraint += cphmd_gen.n_buf * lambdaResidue.lambda.x;
                    }
                    else
                    {
                        lambdaChargeConstraint += lambdaResidue.lambda.x;
                    }
                }
            }
        }
    }

    // iterate until all multiple states and charge constraint are close to desired
    int it = 0;
    while ((it < 10000 && max_deviation > 0.00001)
           || (it < 10000 && abs(lambdaChargeConstraint - lambdaChargeConstraintInitial) > 0.00001))
    {
        // multiple states constraint for all multistate groups
        if (cphmd_gen.haveMultiStateConstraints)
        {
            std::vector<real> multiLambdaGroupState(numberOfMultiStateGroups, 0);
            // loop over all multistate groups
            for (int group = 1; group <= numberOfMultiStateGroups; group++)
            {
                // sum current multistate group lambdas
                int counter = 0;
                for (const auto& lambdaResidue : lambdaResidues)
                {
                    // check if lambda groups belongs to this multiple states group
                    if (lambdaResidue.chargeRestraintGroupNumber == group)
                    {
                        multiLambdaGroupState[group - 1] += lambdaResidue.lambda.x;
                        counter++;
                    }
                }

                const real sum_dt2_per_mass = counter * dt * dt / cphmd_gen.lambda_mass;
                const real multiplier =
                        (multiLambdaGroupState[group - 1] - multiLambdaInitial) / sum_dt2_per_mass;
                multiLambdaGroupState[group - 1] = 0;

                // linear constraint for one group converges in one step
                for (auto& lambdaResidue : lambdaResidues)
                {
                    if (lambdaResidue.chargeRestraintGroupNumber == group)
                    {
                        lambdaResidue.lambda.x -= multiplier * dt * dt / cphmd_gen.lambda_mass;
                        multiLambdaGroupState[group - 1] += lambdaResidue.lambda.x;
                    }
                }
            }

            for (int i = 0; i < numberOfMultiStateGroups; i++)
            {
                deviation[i] = abs(multiLambdaGroupState[i] - multiLambdaInitial);
            }
            max_deviation = *std::max_element(
                    deviation.begin(), deviation.begin() + numberOfMultiStateGroups
                                               - 1); // largest deviation from l_multi_init = 1
        }

        // charge constraint for total sum of lambdas
        if (cphmd_gen.haveChargeConstraints)
        {
            const real sum_dt2_per_mass =
                    (numConstrainedLambdas - 1) * dt * dt / cphmd_gen.lambda_mass
                    + cphmd_gen.n_buf * cphmd_gen.n_buf * dt * dt / cphmd_gen.lambda_mass;
            real lambdaChargeConstraint = 0;

            // sum all lambdas (and collective buffer)
            for (int i = 1; i <= gmx::ssize(lambdaResidues); i++)
            {
                const auto& lambdaResidue = lambdaResidues[i];
                for (int j = 0; j < numConstrainedLambdas; j++)
                {
                    if (i == cphmd_gen.constrainedLambdaIndices[j])
                    {
                        if (lambdaResidue.isBufferResidue)
                        {
                            lambdaChargeConstraint += cphmd_gen.n_buf * lambdaResidue.lambda.x;
                        }
                        else
                        {
                            lambdaChargeConstraint += lambdaResidue.lambda.x;
                        }
                    }
                }
            }

            real multiplier = (lambdaChargeConstraint - lambdaChargeConstraintInitial) / sum_dt2_per_mass;
            lambdaChargeConstraint = 0;

            // linear constraint converges in one step
            for (int i = 1; i <= gmx::ssize(lambdaResidues); i++)
            {
                auto& lambdaResidue = lambdaResidues[i];
                for (int j = 0; j < numConstrainedLambdas; j++)
                {
                    if (i == cphmd_gen.constrainedLambdaIndices[j])
                    {
                        if (lambdaResidue.isBufferResidue)
                        {
                            lambdaResidue.lambda.x -=
                                    cphmd_gen.n_buf * multiplier * dt * dt / cphmd_gen.lambda_mass;
                            lambdaChargeConstraint += cphmd_gen.n_buf * lambdaResidue.lambda.x;
                        }
                        else
                        {
                            lambdaResidue.lambda.x -= multiplier * dt * dt / cphmd_gen.lambda_mass;
                            lambdaChargeConstraint += lambdaResidue.lambda.x;
                        }
                    }
                }
            }
        }
        it++;
    }

    // advance lambda velocities, kinetic energies, temperatures
    // (also for non-constrained groups these now computed again but should not matter)

    for (auto& lambdaResidue : lambdaResidues)
    {
        lambdaResidue.lambda.v = (lambdaResidue.lambda.x - lambdaResidue.lambda.x_old) / dt;
        lambdaResidue.lambda.ekin =
                0.5 * cphmd_gen.lambda_mass * lambdaResidue.lambda.v * lambdaResidue.lambda.v;
        lambdaResidue.lambda.T = 2.0 * (lambdaResidue.lambda.ekin) / BOLTZ;
    }
}


template<int phase>
void do_charge_constraint(const t_inputrec&             ir,
                          const cphmd_general&          cphmd_gen,
                          gmx::ArrayRef<LambdaResidues> ml,
                          const real                    sum_dt2_per_mass)
{
    real       lambda_tot      = 0; // total sum of lambdas
    const real lambda_tot_init = cphmd_gen.lambda_tot_init;
    const real dt              = ir.delta_t;

    // sum current lambda values
    for (const auto& lambdaResidue : ml)
    {
        lambda_tot += lambdaResidue.isBufferResidue ? lambdaResidue.lambda.x * cphmd_gen.n_buf
                                                    : lambdaResidue.lambda.x;
    }

    // linear constraint converges in one step
    const real multiplier = (lambda_tot - lambda_tot_init) / sum_dt2_per_mass;
    for (auto& lambdaResidue : ml)
    {
        lambdaResidue.lambda.x -= lambdaResidue.isBufferResidue
                                          ? cphmd_gen.n_buf * multiplier * dt * dt / cphmd_gen.lambda_mass
                                          : multiplier * dt * dt / cphmd_gen.lambda_mass;
    }

    // phase 0 for v-rescale thermostat
    if (phase == 0)
    {
        for (auto& lambdaResidue : ml)
        {
            lambdaResidue.lambda.v = (lambdaResidue.lambda.x - lambdaResidue.lambda.x_old) / dt;
            lambdaResidue.lambda.ekin =
                    0.5 * cphmd_gen.lambda_mass * lambdaResidue.lambda.v * lambdaResidue.lambda.v;
            lambdaResidue.lambda.T = 2.0 * (lambdaResidue.lambda.ekin) / BOLTZ;
        }
    }

    // phases 1,2 for langevin charge constraint
    if (phase == 1)
    {
        for (auto& lambdaResidue : ml)
        {
            lambdaResidue.lambda.v = (lambdaResidue.lambda.x - lambdaResidue.lambda.x_old_old) / dt;
        }
    }
    if (phase == 2)
    {
        for (auto& lambdaResidue : ml)
        {
            lambdaResidue.lambda.v =
                    lambdaResidue.lambda.v_old_old
                    + 2.0 * (lambdaResidue.lambda.x - lambdaResidue.lambda.x_old_old) / dt;
            lambdaResidue.lambda.T = 2.0 * (lambdaResidue.lambda.ekin) / BOLTZ;
        }
    }
}

/*  Perform T coupling with v-rescale thermostat
 *  - couple all lambdas together
 *  - Nf = number of degrees of freedom = N - 1
 *  - according to Bussi et al JCP (2007) appendix
 *
 *  Returns the change in kinetic energy
 */
real tcouple_vrescale_collective(const cphmd_general&          cphmd_gen,
                                 gmx::ArrayRef<LambdaResidues> ml,
                                 real                          Tref,
                                 real                          tau,
                                 real                          dt,
                                 int64_t                       step,
                                 int64_t                       seed,
                                 bool                          useChargeConstraints,
                                 bool                          useMultiStateConstraits)
{
    gmx::ThreeFry2x64<64> rng(seed, gmx::RandomDomain::Thermostat);

    gmx::NormalDistribution<real> normalDist;

    const real factor = exp(-1.0 * dt / tau);

    // number of degrees of freedom (one less if constraint on)
    int Nf = ml.size();
    if (useChargeConstraints)
    {
        Nf = Nf - 1;
    }
    if (useMultiStateConstraits)
    {
        Nf = Nf - gmx::ssize(cphmd_gen.multiGroupLambdaStates);
    }

    // total kinetic energy
    real Ekin_total = 0;
    for (const auto& lambdaResidue : ml)
    {
        Ekin_total = Ekin_total
                     + half * cphmd_gen.lambda_mass * lambdaResidue.lambda.v * lambdaResidue.lambda.v;
    }

    // reference kinetic energy from desired temperature
    const real Ekin_ref = Nf * 0.5 * Tref * BOLTZ;

    // first random number
    rng.restart(step, 0);
    const real r1 = normalDist(rng);

    // sum of gaussian squared random numbers
    real sum_r2 = 0;
    for (int j = 1; j < Nf; j++)
    {
        const real r = normalDist(rng);
        sum_r2       = sum_r2 + r * r;
    }
    const real Ekin_new = Ekin_total + (1.0 - factor) * (Ekin_ref * (r1 * r1 + sum_r2) / Nf - Ekin_total)
                          + 2.0 * r1 * std::sqrt(Ekin_ref * Ekin_total / Nf * (1.0 - factor) * factor);

    // Analytically Ek_new>=0, but we check for rounding errors (from gromacs coupling.cpp)
    real alpha = 0.0;
    if (Ekin_new > 0)
    {
        alpha = std::sqrt(Ekin_new / Ekin_total);
    }

    // scale kinetic energies and velocities with alpha
    for (auto& lambdaResidue : ml)
    {
        lambdaResidue.lambda.v    = alpha * lambdaResidue.lambda.v;
        lambdaResidue.lambda.x    = lambdaResidue.lambda.x_old + lambdaResidue.lambda.v * dt;
        lambdaResidue.lambda.ekin = alpha * alpha * lambdaResidue.lambda.ekin;
        lambdaResidue.lambda.T    = two * (lambdaResidue.lambda.ekin) / BOLTZ;
    }

    return Ekin_new - Ekin_total;
}

/*
 *  Update lambda and velocity using normal leap-frog
 *  - needs thermostat (collective v-rescale)
 */
void updateLambda(const t_inputrec& ir, LambdaResidues* ml, const cphmd_general& cphmd_gen)
{
    if (!cphmd_gen.isCalibrationRun)
    {
        GMX_RELEASE_ASSERT(cphmd_gen.lambda_mass != 0.0,
                           "Can't have zero mass for non calibration simulations");
        const real inverseMass = 1.0 / cphmd_gen.lambda_mass;

        ml->lambda.v     = ml->lambda.v + inverseMass * (-1.0) * ml->lambda.dvdl * ir.delta_t;
        ml->lambda.ekin  = half * cphmd_gen.lambda_mass * ml->lambda.v * ml->lambda.v;
        ml->lambda.T     = two * (ml->lambda.ekin) / BOLTZ;
        ml->lambda.x_old = ml->lambda.x;
        ml->lambda.x     = ml->lambda.x + ml->lambda.v * ir.delta_t;
    }
}

} // namespace

ConstantPH::ConstantPH(const t_inputrec& ir, int natoms, const std::string& inputFileName) :
    inputFileName_(inputFileName),
    eLambdaThermostat_(ir.eLambdaThermostat),
    useChargeConstraints_(ir.bLambdaChargeConstraints),
    useMultiStateConstraits_(ir.bLambdaMultiStateConstraints),
    isCalibrationRun_(ir.bLambdaIsCalibration)
{
    GMX_RELEASE_ASSERT(ir.lambda_dynamics, "We should only set up ConstantPH with lambda dynamics");

    ml_ = std::make_unique<LambdaResidueCollection>();

    cphmd_gen_ = std::make_unique<cphmd_general>();
    init_constantph(ml_.get(), cphmd_gen_.get(), inputFileName_, useChargeConstraints_,
                    useMultiStateConstraits_, isCalibrationRun_);

    potential_.resize(natoms);

    // copy initial charges of all atoms
    // TODO: This does not work with domain decomposition

    // Create indices
    for (const auto& lambdaResidue : ml_->residues_)
    {
        int i = 0;
        for (const auto& atom : lambdaResidue.atoms)
        {
            lambdaAtoms_.push_back(atom);
            lambdaAtomsIndex_.push_back(atom);
            i++;
        }
    }
    std::vector<int> unsortedLambdaAtoms(lambdaAtoms_);
    std::sort(lambdaAtoms_.begin(), lambdaAtoms_.end());
    std::vector<int> sortedLambdaAtoms;
    for (size_t i = 0; i < lambdaAtomsIndex_.size(); ++i)
    {
        auto it = find(unsortedLambdaAtoms.begin(), unsortedLambdaAtoms.end(), lambdaAtoms_[i]);
        int  atomIndex = std::distance(unsortedLambdaAtoms.begin(), it);
        sortedLambdaAtoms.push_back(lambdaAtomsIndex_[atomIndex]);
    }
    lambdaAtomsIndex_ = sortedLambdaAtoms;
}

ConstantPH::~ConstantPH()
{
    // TODO: We should actually clean up or use only C++ constructs
}

void ConstantPH::setLambdaCharges(t_mdatoms* mdatoms)
{
    // interpolate the charges for atoms that are part of lambda group
    for (const auto& lambdaResidue : ml_->residues_)
    {
        int atomInResidue = 0;
        for (const auto& atom : lambdaResidue.atoms)
        {
            mdatoms->chargeA[atom] = (1 - lambdaResidue.lambda.x) * lambdaResidue.chargeA[atomInResidue]
                                     + lambdaResidue.lambda.x * lambdaResidue.chargeB[atomInResidue];
            atomInResidue += 1;
        }
    }
}

/* TODO
 * loop over lambda groups
 *     lambda = do_lambdadyn(...)
 * output lambda value, temperature etc.
 *
 * Returns the change in kinetic energy due to T-coupling
 */
real ConstantPH::updateLambdas(const t_inputrec& ir, const double t, const int64_t step)
{
    // dvdl for each lambda group -- now use this in force calculation
    std::vector<real> dvdl(ml_->residues_.size(), 0);
    {
        int residueIndex = 0;
        int atomIndex    = 0;
        for (const auto& lambdaResidue : ml_->residues_)
        {
            int atomInResidue = 0;
            for (const gmx_unused auto& atom : lambdaResidue.atoms)
            {
                dvdl[residueIndex] +=
                        potential_[lambdaAtomsIndex_[atomIndex]]
                        * (lambdaResidue.chargeB[atomInResidue] - lambdaResidue.chargeA[atomInResidue]);
                atomInResidue++;
                atomIndex++;
            }
            residueIndex++;
        }
    }
    // for multistates put lambdas to table
    std::vector<real> lambdaValues;
    if (useMultiStateConstraits_)
    {
        for (const auto& lambdaResidue : ml_->residues_)
        {
            lambdaValues.emplace_back(lambdaResidue.lambda.x);
        }
    }
    {
        int i = 0;
        for (auto& lambdaResidue : ml_->residues_)
        // compute forces acting on each lambda
        {
            compute_forces(*cphmd_gen_, &lambdaResidue, ml_->residues_, dvdl[i], i);
            i++;
        }
    }

    /* Thermostat is controlled by mdp option.
     */

    // if we have thermostat = langevin
    const real gamma     = cphmd_gen_->tau_lambda;
    real       deltaEkin = 0;
    switch (eLambdaThermostat_)
    {
        case (eLambdaTcLangevin):
        {
            // lambda update if we have charge constraint on (two phases and constraints)
            if (useChargeConstraints_)
            {
                for (auto& lambdaResidue : ml_->residues_)
                {
                    updateLambdaLD<1>(ir, *cphmd_gen_, &lambdaResidue, gamma);
                }

                const real sum_dt2_per_mass =
                        ((ml_->residues_.size() - 1) * (ir.delta_t) * (ir.delta_t)) / cphmd_gen_->lambda_mass
                        + (cphmd_gen_->n_buf * cphmd_gen_->n_buf * ir.delta_t * ir.delta_t)
                                  / cphmd_gen_->m_buf;
                do_charge_constraint<1>(ir, *cphmd_gen_, ml_->residues_, sum_dt2_per_mass);

                for (auto& lambdaResidue : ml_->residues_)
                {
                    updateLambdaLD<2>(ir, *cphmd_gen_, &lambdaResidue, gamma);
                }

                do_charge_constraint<2>(ir, *cphmd_gen_, ml_->residues_, sum_dt2_per_mass);
            }

            // if charge constraint is not applied, normal Langevin update
            else
            {
                for (auto& lambdaResidue : ml_->residues_)
                {
                    updateLambdaLD<3>(ir, *cphmd_gen_, &lambdaResidue, gamma);
                }
            }
        }
        break;
        case (eLambdaTcVRESCALE):

        {
            for (auto& lambdaResidue : ml_->residues_)
            {
                updateLambda(ir, &lambdaResidue, *cphmd_gen_);
            }

            if (cphmd_gen_->lambda_mass != 0)
            {
                deltaEkin = tcouple_vrescale_collective(
                        *cphmd_gen_, ml_->residues_, cphmd_gen_->reference_temperature,
                        cphmd_gen_->tau_lambda, ir.delta_t, step, ir.ld_seed, useChargeConstraints_,
                        useMultiStateConstraits_);
            }

            // multiple states
            if (useMultiStateConstraits_ || useChargeConstraints_)
            {
                do_constraints(*cphmd_gen_, ml_->residues_, ir.delta_t);
            }
        }
        break;
        default: GMX_RELEASE_ASSERT(false, "How did we end up here?");
    }

    // print results for one time step
    for (const auto& lambdaResidue : ml_->residues_)
    {
        // only print out every nst_lambda steps
        int step = std::floor(t / (&ir)->delta_t);
        if (step % cphmd_gen_->nst_lambda == 0)
        {
            fprintf(lambdaResidue.out, "%f %g %g %g %g %g %g %g %g\n ", t, lambdaResidue.lambda.x,
                    lambdaResidue.lambda.dvdl, lambdaResidue.lambda.T, lambdaResidue.lambda.v,
                    lambdaResidue.lambda.dvdl_pot, lambdaResidue.lambda.dvdl_ref,
                    lambdaResidue.lambda.dvdl_dwp, lambdaResidue.lambda.dvdl_ph);
        }
    }

    return deltaEkin;
}
