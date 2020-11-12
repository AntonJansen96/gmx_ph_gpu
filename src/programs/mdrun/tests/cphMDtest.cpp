/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * This implements basic constant pH MD sanity tests.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief A basic cpHMD runner.
 */
class CphMDTest : public MdrunTestFixture, public ::testing::WithParamInterface<int>
{
public:
    //! Runs the test with the given inputs
    void runTest();
};

void CphMDTest::runTest()
{
    runner_.ndxFileName_ = "";
    ASSERT_EQ(0, runner_.callGrompp());

    CommandLine commandLine;
    std::string cphMDInputFile = TestFileManager::getInputFilePath("GLU_constant_ph_input.dat");
    commandLine.append("mdrun");
    commandLine.addOption("-cpHMD", cphMDInputFile.c_str());
    ASSERT_EQ(0, runner_.callMdrun(commandLine));
}

TEST_F(CphMDTest, RunsGLUTestCase)
{
    runner_.useTopGroAndNdxFromDatabase("GLU");
    const std::string mdpFileContents = formatString(R"(
; Run control 
integrator               = md                ; md, md-vv, steep, cg ...
tinit                    = 0                 ; starting time for run (0)
dt                       = 0.002             ; timestep [ps] (0.001)
nsteps                   = 200             ; number of steps (0)
init-step                = 0                 ; starting step (0)
comm-mode                = Linear            ; mode for center of mass motion removal
nstcomm                  = 10                ; number of steps for center of mass motion  removal (100)

; Output control options 

nstxout                  = 100               ; output frequency for coords x (0)
nstvout                  = 100               ; output frequency for velocities v (0)
nstfout                  = 0                 ; output frequency for forces f (0)
nstlog                   = 100               ; output frequency for energies to log file (1000)
nstcalcenergy            = -1                ; frequency for calculating energies (100)
nstenergy                = 100               ; output frequency for energies to energy file (1000)
nstxout-compressed       = 0                 ; output frequency for writing x using lossy compression (0)

; Neighbor searching parameters
cutoff-scheme            = Verlet            ; verlet or group
nstlist                  = 10               ; neighbor list update frequency (10)
ns-type                  = grid              ; ns algorithm (simple or grid)
pbc                      = xyz               ; periodic boundary conditions: xyz, no, xy
periodic_molecules       = no                ; no = finite molecules, yes = molecules couple through pbc
verlet-buffer-tolerance  = 0.005             ; only with verlet, sets indirectly rlists (0.005)
rlist                    = 1.0               ; nblist cut-off (ignored if verlet, buffer sets it) [nm] (1) 



; Options for electrostatics 

coulombtype              = PME               ; method for doing electrostatics (PME standard)
coulomb-modifier         = potential-shift-verlet ; use with verlet, selects none for group scheme
rcoulomb-switch          = 0                 ; where to start switching coulomb potential if it's on [nm] (0)
rcoulomb                 = 1.0               ; distance for coulomb cut-off [nm] (1)

; Options for van der Waals 

vdwtype                  = cut-off           ; method for doing van der Waals  
rvdw                     = 1.0               ; distance for LJ cut-off [nm] (1)

; Temperature coupling

tcoupl                   = berendsen         ; temperature coupling (no, berendsen, v-rescale, nose-hoover)
nsttcouple               = -1                ; frequency for coupling temperature (-1, sets equal to nstlist)
tc-grps                  = System            ; groups to couple separately (water, protein)
tau-t                    = 0.1               ; time constant [ps] for each group
ref-t                    = 300               ; reference temperature [K] for each group   



; Pressure coupling
 
pcoupl                   = berendsen         ; pressure coupling (no, berendse, parrinello-rahman) 
Pcoupltype               = isotropic         ; usually isotropic, semiisotropic
nstpcouple               = -1                ; frequency for coupling pressure (-1, sets equal to nstlist)
tau-p                    = 1.0               ; time constant [ps]
compressibility          = 4.5e-5            ; compressibility [1/bar]
ref-p                    = 1.0               ; reference P [bar]
refcoord_scaling         = No                ; scaling of reference coordinates, No, All or COM
; Generate velocities 

gen-vel                  = yes               ; if changed, also change 'continuation'
gen-temp                 = 300               ; temperature for Maxwell distribution [K] (300)

; Options for bonds   
 
constraints              = h-bonds           ; none, h-bonds, all-bonds
constraint-algorithm     = lincs             ; type of constraint algorithm (lincs, shake)

; Free energy control 

free-energy              = no              ; yes = interpolated between topology A and B
init-lambda              = 0                 ; starting value for lambda, only for slow growth
delta-lambda             = 0                 ; increment per time step for lambda
init-lambda-state        = -1                ; which column of lambda vector should be used (-1)

; lambda values for which dh values are determined, values between 0 and 1 (fep for both coul and vdw)
fep-lambdas              =                   
coul-lambdas             = 
vdw-lambdas              =
calc-lambda-neighbors    = 1                 ; number of neighbor lambdas for which dh will be calculated (1)
sc-alpha                 = 0                 ; soft-core alpha parameter, (0) results in linear interpolation
sc-r-power               = 6                 ; power of radial term in soft-core, values 6 or 48 (6)
sc-coul                  = no                ; whether to apply soft-core to coulomb (no)
sc-power                 = 0                 ; power for lambda in soft-core function, values 1 and 2 (0)
sc-sigma                 = 0.3               ; soft core sigma [nm] (0.3)
couple-moltype           =                   ; molecule type for solvation of coupling energies
couple-lambda0           = vdw-q             ; which interactions are on at lambda=0 (vdw-q/vdw/q/none)
couple-lambda1           = vdw-q             ; which interactions are on at lambda=1 (vdw-q/vdw/q/none)
couple-intramol          = no                ; replace intra-molecular interactions by exclusions
nstdhdl                  = 10                ; frequency for writing dH/dlambda to dhdl.xvg (100)
dhdl-derivatives         = yes               ; write out derivatives of hamiltonian
dhdl-print-energy        = no                ; include potential or total energy in dhdl file
separate-dhdl-file       = yes               ; write dhdls to a separate file, by default dhdl.xvg
dh_hist_size             = 0                 ; nonzero value speifies size of histogram for dhdl values
dh_hist_spacing          = 0.1               ; bin width of the histograms in energy units

lambda-dynamics          = yes
    )");

    runner_.useStringAsMdpFile(mdpFileContents);
    runTest();
}

} // namespace
} // namespace test
} // namespace gmx
