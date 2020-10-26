/*
 * constant_ph.h
 *
 */

#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/random/threefry.h"
#include <stdio.h>
#include <config.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
using namespace std;

/* Structure for one lambda group */
struct t_lambdarec {
    /* x is for lambda coordinate */
    real  x0;                     /* initial lambda */
    real  x;                      /* lambda coordinate */
    real  x_old;                  /* lambda of previous step */
    real  x_old_old;              /* lambda of previous step */
    real  v;                      /* velocity of lambda */
    real  v_old;                  /* velocity of previous step */
    real  v_old_old;              /* lambda of previous step */
    real  ekin;                   /* kinetic energy of lambda */
    real  T;                      /* temperature */
    int   tau;                    /* tcouple time constant */
    real  m;                      /* mass */
    real  bar;                    /* double well barrier */
    real  dvdl;                   /* dV/dlambda, force acting on lambda */
    real  dvdl_old;               /* dV/dlambda, force acting on lambda */
    real *lambda_dwp;             /* double well potential array*/
    real dvdl_pot;  /* dvdl from environment */
    real dvdl_ph;   /* dvdl from ph */
    real dvdl_dwp;  /* dvdl from double well potential */
    real dvdl_ref;  /* dvdl from reference free energy */
};

/* Structure for all lambda groups */
struct multi_lambda{
    FILE                 *out;             /* output file for lambda dynamics data */
    struct t_lambdarec   *lambda;          /* lambda structure */
    real                 *chargeA;         /* charge A */
    real                 *chargeB;         /* charge B */
    std::string           residue_name;    /* name of residue */
    int                   residue_index;   /* reference to residues in cphmd_general*/
    int                   group_number;    /* group number for multiple states */
    int                   n_atoms;         /* number of atoms in lambda group */
    int                  *atoms;           /* indexes of atoms of lambda group */
    struct multi_lambda  *next;            /* pointer to next lambda group */
};

/* Structure for general and residue specific data */
struct cphmd_general{
    real                     ph_value;          /* value of pH */
    std::string              charge_constraint; /* takes values "no", "yes"*/
    std::string              multistate_constraint;  /* takes values "no", "yes" */ 
    int                      n_multigroups;     /* number of different constraint groups*/
	int                      *n_states;         /* array for how many states in each multigroup lambda */
    int                      n_buf;             /* number of (collective) buffers */
    real                     m_buf;             /* mass of buffer*/
	int                      n_constrained_lambdas; /* how many lmabdas we have for charge constraint */
	int                     *constrained_lambdas; /* indexes for lambdas one wants to constrain */
    real                     lambda_tot_init;   /* sum of lambdas for charge constraint */
    real                     T_lambda;        /* T of lambda */
    real                     m_lambda;        /* mass of lambda*/
    int                      nst_lambda;      /* output frequency*/
    real                     tau;             /* time constant / langevin friction rate */
    std::string              thermo;          /* 'langevin' or 'v-rescale' choose how to do it*/
    int                      nr_res;          /* number of residues */
    int                      nr_lg;           /* number of lambda groups */
    int                      nr_lambda_atoms; /* total number of lambda atoms*/
    std::vector<std::string> residue_names;   /* array for residue names */
    real                    *pKa_values;      /* array for reference pKa values */
    int                     *n_coeffs;        /* array for n_coeffs */
    real                    *dvdl_coeffs;     /* array (1D) for dvdl_coeffs (now max 10 coeffs per fit) */
	real                    *dvdl_coeffs_other;     /* another array for dvdl_coeffs for multiple state edges */
};
extern void updateLambdaLD(const t_inputrec *ir, struct multi_lambda *ml, struct cphmd_general *cphmd_gen, real T_ref, real gamma, int phase);
extern void updateLambda(const t_inputrec *ir, struct multi_lambda *ml, struct cphmd_general *cphmd_gen);
extern void compute_forces(struct cphmd_general *cphmd_gen, struct multi_lambda *ml, real dvdl, real lambdas[], int l_index);
extern void init_constantph(struct multi_lambda *ml, struct cphmd_general *cphmd_gen);
extern void do_constraints(struct cphmd_general *cphmd_gen, struct multi_lambda *ml, real dt);
extern void do_charge_constraint(struct cphmd_general *cphmd_gen, struct multi_lambda *ml, real sum_dt2_per_mass, real dt, int phase);
extern void do_multiple_states_constraint(struct cphmd_general *cphmd_gen, struct multi_lambda *ml, real dt, int group);
extern struct multi_lambda  *ml;
extern struct cphmd_general *cphmd_gen;
