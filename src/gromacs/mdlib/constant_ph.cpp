/*
 * Constant pH MD to run with Gromacs
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

real tcouple_vrescale_collective(struct cphmd_general *cphmd_gen, struct multi_lambda *ml, real Tref, real tau, real dt, int64_t step, int64_t seed);


ConstantPH::ConstantPH(const t_inputrec &ir,
                       const t_mdatoms  &mdatoms)
{
    GMX_RELEASE_ASSERT(ir.lambda_dynamics, "We should only set up ConstantPH with lambda dynamics");

    ml_           = new multi_lambda;

    ml_temp_      = ml_; // save first lambda in ml for later use

    cphmd_gen_    = std::make_unique<cphmd_general>();
    init_constantph(ml_, cphmd_gen_.get());

    ml_ = ml_temp_; // go back to beginning of ml (first lambda in ml)

    potential_.resize(mdatoms.nr);

    // copy initial charges of all atoms
    // TODO: This does not work with domain decomposition
    chargeA_.resize(mdatoms.nr);
    chargeB_.resize(mdatoms.nr);
    for (int j = 0; j < mdatoms.nr; j++)
    {
        chargeA_[j] = mdatoms.chargeA[j];
        chargeB_[j] = mdatoms.chargeA[j]; // same charges
    }

    // Create indices
    for (int i = 0; i < cphmd_gen_->nr_lg; i++)
    {
        for (int j = 0; j < ml_->n_atoms; j++)
        {
            lambdaAtoms_.push_back(ml_->atoms[j]);
            lambdaAtomsIndex_.push_back(i);
            chargeA_[ml_->atoms[j]] = ml_->chargeA[j];
            chargeB_[ml_->atoms[j]] = ml_->chargeB[j];
        }
        ml_ = ml_->next;
    }
	ml_ = ml_temp_; // go back to beginning of ml (first lambda in ml)
    dvdl_.resize(cphmd_gen_->nr_lg);

	// sort lambdaAtoms and lambdaAtomsIndex
	std::vector<int> lA_notsorted;
	lA_notsorted = lambdaAtoms_;
    std::vector<int> lAInd_sorted;
	int ind;

	sort(lambdaAtoms_.begin(), lambdaAtoms_.end());

    for (std::size_t i = 0; i < lambdaAtomsIndex_.size(); ++i)
    {
		auto it = find(lA_notsorted.begin(),lA_notsorted.end(),lambdaAtoms_[i]);
		ind = std::distance(lA_notsorted.begin(), it);
		lAInd_sorted.push_back(lambdaAtomsIndex_[ind]);
    }
	lambdaAtomsIndex_ = lAInd_sorted;

	for (std::size_t i = 0; i < lambdaAtomsIndex_.size(); ++i)
	{
        fprintf(stderr, "atom %i lambda group %i chargeA %g chargeB %g  \n",lambdaAtoms_[i],lambdaAtomsIndex_[i]+1,chargeA_[lambdaAtoms_[i]],chargeB_[lambdaAtoms_[i]]);
	}
}

ConstantPH::~ConstantPH()
{
    // TODO: We should actually clean up or use only C++ constructs
}

void ConstantPH::setLambdaCharges(t_mdatoms *mdatoms)
{
    ml_ = ml_temp_; // return to first lambda group

    // interpolate the charges for atoms that are part of lambda group

    //int pos = 0;
    //for (int i = 0; i < cphmd_gen_->nr_lg; i++)
    //{
	//    for (int j = 0; j < ml_->n_atoms; j++)
    //    {
    //        const int atomIndex = lambdaAtoms_[pos];
    //        mdatoms->chargeA[atomIndex] =
    //            (1 - ml_->lambda->x)*chargeA_[atomIndex] + ml_->lambda->x*chargeB_[atomIndex];
    //        pos += 1;
    //    }
    //    ml_ = ml_->next;
    //}
    //ml_ = ml_temp_;

    // does it matter if we don't use lambdaAtoms ?
    for (int i = 0; i < cphmd_gen_->nr_lg; i++)
    {
        for (int j = 0; j < ml_->n_atoms; j++)
        {
            const int atomIndex = ml_->atoms[j];
            mdatoms->chargeA[atomIndex] =
                (1 - ml_->lambda->x)*chargeA_[atomIndex] + ml_->lambda->x*chargeB_[atomIndex];
        }
        ml_ = ml_->next;
    }
    ml_ = ml_temp_;

}

FILE *out_T;

/* TODO
 * loop over lambda groups
 *     lambda = do_lambdadyn(...)
 * output lambda value, temperature etc.
 *
 * Returns the change in kinetic energy due to T-coupling
 */
real ConstantPH::updateLambdas(const t_inputrec &ir,
                               const double      t,
							   const int64_t     step)
{
    ml_ = ml_temp_;

    // dvdl for each lambda group -- now use this in force calculation
    std::fill(dvdl_.begin(), dvdl_.end(), 0);
    for (std::size_t i = 0; i < lambdaAtoms_.size(); ++i)
    {
       dvdl_[lambdaAtomsIndex_[i]] += potential_[lambdaAtoms_[i]]*(chargeB_[lambdaAtoms_[i]]-chargeA_[lambdaAtoms_[i]]);
	   //fprintf(stderr,"%zu %i %i %g\n",i,lambdaAtomsIndex_[i],lambdaAtoms_[i],dvdl_[lambdaAtomsIndex_[i]]);   
    }

    // for multistates put lambdas to table
    real lambdas [cphmd_gen_->nr_lg];
	if(cphmd_gen_.get()->multistate_constraint.compare("yes")==0)
	{
        for (int i = 0; i < cphmd_gen_->nr_lg; i++)
        {
            lambdas[i] = ml_->lambda->x;
            ml_ = ml_->next;
        }
	    ml_ = ml_temp_;
    }

    // compute forces acting on each lambda
    for (int i = 0; i < cphmd_gen_->nr_lg; i++)
    {
        compute_forces(cphmd_gen_.get(), ml_, dvdl_[i], lambdas, i);
		//fprintf(stderr,"dvdl %i %g\n",i, dvdl_[i]);
        ml_ = ml_->next;
    }
	ml_ = ml_temp_;
			
    real sum_dt2_per_mass;

    /* in constant_ph_input.dat we have "thermostat" option:
     * thermostat = langevin
     * thermostat = v-rescale
     */
	
    // if we have thermostat = langevin
	real gamma = cphmd_gen_->tau;
	int phase;
    if(cphmd_gen_.get()->thermo.compare("langevin")==0)
    {
	    // lambda update if we have charge constraint on (two phases and constraints)
	    if(cphmd_gen_.get()->charge_constraint.compare("yes")==0)
	    {
		    phase = 1;    // phase = 1 update, no noise and friction yet
	        for (int i = 0; i < cphmd_gen_->nr_lg; i++)
	        {
		        updateLambdaLD(&ir, ml_, cphmd_gen_.get(), cphmd_gen_->T_lambda, gamma, phase);
		        ml_ = ml_->next;
	        }
            ml_ = ml_temp_;

            sum_dt2_per_mass = (cphmd_gen_->nr_lg-1)*(&ir)->delta_t*(&ir)->delta_t/cphmd_gen_->m_lambda + cphmd_gen_->n_buf*cphmd_gen_->n_buf*(&ir)->delta_t*(&ir)->delta_t/cphmd_gen_->m_buf;
            do_charge_constraint(cphmd_gen_.get(), ml_, sum_dt2_per_mass, (&ir)->delta_t, phase);  
	        ml_ = ml_temp_;

	        phase = 2;    // phase 2: add langevin friction and noise
	        for (int i = 0; i < cphmd_gen_->nr_lg; i++)
	        {
	        	updateLambdaLD(&ir, ml_, cphmd_gen_.get(), cphmd_gen_->T_lambda, gamma, phase);
	        	ml_ = ml_->next;
	        }
	        ml_ = ml_temp_;

    	    do_charge_constraint(cphmd_gen_.get(), ml_, sum_dt2_per_mass, (&ir)->delta_t, phase);
    	    ml_ = ml_temp_;
	    }
    
        // if charge constraint is not applied, normal Langevin update
	    if(cphmd_gen_.get()->charge_constraint.compare("off")==0)
	    {
    	    phase = 3;
    	    for (int i = 0; i < cphmd_gen_->nr_lg; i++)
    	    {
    	    	updateLambdaLD(&ir, ml_, cphmd_gen_.get(), cphmd_gen_->T_lambda, gamma, phase);
    	    	ml_ = ml_->next;
    	    }
            ml_ = ml_temp_;
        }
    }

    real deltaEkin = 0;

    // if we have thermostat = v-rescale
    if(cphmd_gen_.get()->thermo.compare("v-rescale")==0)
    {
        for (int i = 0; i < cphmd_gen_->nr_lg; i++)
    	{
    		updateLambda(&ir, ml_, cphmd_gen_.get());
    		ml_ = ml_->next;
    	}
    	ml_ = ml_temp_;
        
        if (cphmd_gen_->m_lambda != 0)
    	{
            deltaEkin =
                tcouple_vrescale_collective(cphmd_gen_.get(), ml_, cphmd_gen_->T_lambda, cphmd_gen_->tau, (&ir)->delta_t, step, (&ir)->ld_seed);
			
			// general function to loop over constraints
	    	if(cphmd_gen_.get()->multistate_constraint.compare("yes")==0 || cphmd_gen_.get()->charge_constraint.compare("yes")==0)
	    	{
				do_constraints(cphmd_gen_.get(), ml_, (&ir)->delta_t);
			}
    	}
    	ml_ = ml_temp_;
    }

    // compute the average temperature for lambdas (less degrees of freedom if any constraint applied)
    real Ekin_tot = 0;
    real T = 0;
    for (int i = 0; i < cphmd_gen_->nr_lg; i++)
    {
    	Ekin_tot = Ekin_tot + ml_->lambda->ekin;
    	ml_ = ml_->next;
    }

    int Nf = cphmd_gen_->nr_lg;
    if(cphmd_gen_->charge_constraint.compare("yes")==0)
    {
    	Nf = Nf - 1;
    }
    if(cphmd_gen_->multistate_constraint.compare("yes")==0)
    {
    	Nf = Nf - cphmd_gen_->n_multigroups;
    }
	T = 2*Ekin_tot/(Nf*BOLTZ);
    ml_ = ml_temp_;

    if(step % cphmd_gen_->nst_lambda == 0)
    {
    	fprintf(out_T, "%f %g \n", t, T);
    }

    // print results for one time step
    for (int i = 0; i < cphmd_gen_->nr_lg; i++)
    {
        // only print out every nst_lambda steps
        int step = std::floor(t/(&ir)->delta_t);
        if(step % cphmd_gen_->nst_lambda == 0)
        {
			fprintf(stderr, "time %.6f lambda %.4f dvdl %.4f T %.4f\n", t, ml_->lambda->x, ml_->lambda->dvdl,ml_->lambda->T);

            fprintf(ml_->out, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n ", t, ml_->lambda->x, ml_->lambda->dvdl, ml_->lambda->T, ml_->lambda->v, ml_->lambda->dvdl_pot, ml_->lambda->dvdl_ref,  ml_->lambda->dvdl_dwp,  ml_->lambda->dvdl_ph);
        }
        ml_ = ml_->next;
    }
    ml_ = ml_temp_;

    return deltaEkin;
}

using namespace std;

real gaussdist(gmx::DefaultRandomEngine *rng, real sigma);
void init_lambda_dwp(real *lambda_dwp, real barrier);

int MAX_N_DVDL_COEFFS = 10;

/* For the velocity initialization */
real gaussdist(gmx::DefaultRandomEngine *rng, real sigma)
{
    real r = 2.0, x, y;
    gmx::UniformRealDistribution<real> uniformDist;

    do
    {
        x = 2.0 * uniformDist(*rng) - 1.0;
        y = 2.0 * uniformDist(*rng) - 1.0;
        r = x * x + y * y;
    }
    while (r > 1.0 || r == 0.0);
    r = x * sqrt(-2.0 * std::log(r) / r);
    r = r * sigma;
    return (r);
}

/*  TODO: Change to standard GROMACS input  */
void init_constantph(struct multi_lambda *ml, struct cphmd_general *cphmd_gen)
{
    fprintf(stderr, "\nConstant pH MD initialization: \n\n");

    ifstream     myfile( "constant_ph_input.dat" );
    string       temp[3];
    string       line;
    stringstream ss;

    getline( myfile, line );
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->ph_value;
    cout << temp[0] <<" "<< cphmd_gen->ph_value << "\n";
    ss.clear();

    getline( myfile, line );
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->nr_res;
    cout << temp[0] <<" "<< cphmd_gen->nr_res << "\n";
    ss.clear();

    getline( myfile, line );
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->nr_lg;
    cout << temp[0] <<" "<< cphmd_gen->nr_lg << "\n";
    ss.clear();

    getline( myfile, line );

    getline( myfile, line );
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->m_lambda;
    cout << temp[0] <<" "<< cphmd_gen->m_lambda << "\n";
    ss.clear();

    getline( myfile, line );
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->T_lambda;
    cout << temp[0] <<" "<< cphmd_gen->T_lambda << "\n";
    ss.clear();

    getline( myfile, line );
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->tau;
    cout << temp[0] <<" "<< cphmd_gen->tau << "\n";
    ss.clear();

    getline( myfile, line );
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->thermo;
    cout << temp[0] <<" "<< cphmd_gen->thermo << "\n";
    ss.clear();

    getline( myfile, line );
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->nst_lambda;
    cout << temp[0] <<" "<< cphmd_gen->nst_lambda << "\n";
    ss.clear();
	
    getline( myfile, line );
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->multistate_constraint;
    cout << temp[0] <<" "<< cphmd_gen->multistate_constraint << "\n";
    ss.clear();
	
    getline( myfile, line );
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->n_multigroups;
    cout << temp[0] <<" "<< cphmd_gen->n_multigroups << "\n";
    ss.clear();
	
	cphmd_gen->n_states = new int[cphmd_gen->n_multigroups];
	
	int sum_states = 0;
	
    getline( myfile, line );
    ss << line;
    ss >> temp[0] >> temp[1];
    cout << temp[0] << " " << temp[1] << " ";
    for (int k = 0; k < cphmd_gen->n_multigroups; k++)
    {
        ss >> cphmd_gen->n_states[k];
        cout << cphmd_gen->n_states[k] << " ";
		sum_states += cphmd_gen->n_states[k];
    }
    ss.clear();
	cout << "\n";
	
    getline( myfile, line );
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->charge_constraint;
    cout << temp[0] <<" "<< cphmd_gen->charge_constraint << "\n";
    ss.clear();

    getline( myfile, line );
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->n_buf;
    cout << temp[0] <<" "<< cphmd_gen->n_buf << "\n";
    ss.clear();

    getline( myfile, line );
    ss << line;
    ss >> temp[0] >> temp[1] >> cphmd_gen->m_buf;
    cout << temp[0] <<" "<< cphmd_gen->m_buf << "\n";
    ss.clear();

    if (cphmd_gen->charge_constraint.compare("yes")==0) 
	{
		cphmd_gen->n_constrained_lambdas = cphmd_gen->nr_lg + cphmd_gen->n_multigroups - sum_states;
	    cphmd_gen->constrained_lambdas = new int[cphmd_gen->n_constrained_lambdas];
    }
	else cphmd_gen->n_constrained_lambdas = 0;
	
	
    getline( myfile, line );
    ss << line;
    ss >> temp[0] >> temp[1];
    cout << temp[0] << " " << temp[1] << " ";
    for (int k = 0; k < cphmd_gen->n_constrained_lambdas; k++)
    {
        ss >> cphmd_gen->constrained_lambdas[k];
        cout << cphmd_gen->constrained_lambdas[k] << " ";
    }
    ss.clear();
	
	cout << "\n\n";

    cphmd_gen->residue_names.resize(cphmd_gen->nr_res);
    cphmd_gen->pKa_values  = new real[cphmd_gen->nr_res];
    cphmd_gen->n_coeffs    = new int[cphmd_gen->nr_res];

    cphmd_gen->dvdl_coeffs = new real[cphmd_gen->nr_res*MAX_N_DVDL_COEFFS];
	cphmd_gen->dvdl_coeffs_other = new real[cphmd_gen->nr_res*MAX_N_DVDL_COEFFS];

    getline( myfile, line ); // empty line
	
    /* Loop over residue inputs */
    for (int j = 0; j < cphmd_gen->nr_res; j++)
    {
        getline( myfile, line );
        ss << line;
        ss >> temp[0] >> temp[1] >> cphmd_gen->residue_names[j];
        cout << temp[0] <<" "<< cphmd_gen->residue_names[j] << "\n";
        ss.clear();
		
        getline( myfile, line );
        ss << line;
        ss >> temp[0] >> temp[1] >> cphmd_gen->pKa_values[j];
        /* Set reference pKa equal to pH */
        if (cphmd_gen->pKa_values[j]==0)
        {
            cphmd_gen->pKa_values[j] = cphmd_gen->ph_value;
        }
        cout << temp[0] <<" "<< cphmd_gen->pKa_values[j] << "\n";
        ss.clear();

        getline( myfile, line );
        ss << line;
        ss >> temp[0] >> temp[1] >> cphmd_gen->n_coeffs[j];
        cout << temp[0] <<" "<< cphmd_gen->n_coeffs[j] << "\n";
        ss.clear();

        getline( myfile, line );
        ss << line;
        ss >> temp[0] >> temp[1];
        cout << temp[0] << " " << temp[1] << " ";
        for (int k = 0; k < cphmd_gen->n_coeffs[j]; k++)
        {
            ss >> cphmd_gen->dvdl_coeffs[j*MAX_N_DVDL_COEFFS+k];
            cout << cphmd_gen->dvdl_coeffs[j*MAX_N_DVDL_COEFFS+k] << " ";
        }
        ss.clear();
		cout << "\n";
		
		// if we need another dvdl array for three state multistate with edges ref dvdl
		getline( myfile, line );
		ss << line;
		if (line.empty()) 
		{
			cout << "\n";
			ss.clear();
		}
		else
		{
	        ss >> temp[0] >> temp[1];
	        cout << temp[0] << " " << temp[1] << " ";
	        for (int k = 0; k < cphmd_gen->n_coeffs[j]; k++)
	        {
	            ss >> cphmd_gen->dvdl_coeffs_other[j*MAX_N_DVDL_COEFFS+k];
	            cout << cphmd_gen->dvdl_coeffs_other[j*MAX_N_DVDL_COEFFS+k] << " ";
	        }
			ss.clear();
			cout << "\n\n";
			getline( myfile, line ); // empty line
		}
    }

    /* Loop over lambda groups ml->next */
    int pos = 0;
    cphmd_gen->lambda_tot_init = 0;

    for (int j = 0; j < cphmd_gen->nr_lg; j++)
    {
        getline( myfile, line );
        ss << line;
        ss >> temp[0] >> temp[1] >> ml->residue_name;
        cout << temp[0] <<" "<< ml->residue_name << "\n";
        ss.clear();

        // check which residue in cphmd_gen the lambda group relates to
        for (int k = 0; k < cphmd_gen->nr_res; k++)
        {
            if (ml->residue_name.compare(cphmd_gen->residue_names[k]) == 0)
            {
                ml->residue_index = k;
                break;
            }
        }

        getline( myfile, line );
        ss << line;
        ss >> temp[0] >> temp[1] >> ml->group_number;
        cout << temp[0] <<" "<< ml->group_number << "\n";
        ss.clear();

        struct t_lambdarec *lrec = new t_lambdarec;
        ml->lambda = lrec;

        getline( myfile, line );
        ss << line;
        ss >> temp[0] >> temp[1] >> ml->lambda->x0;
        cout << temp[0] <<" "<< ml->lambda->x0 << "\n";
        ss.clear();

        getline( myfile, line );
        ss << line;
        ss >> temp[0] >> temp[1] >> ml->lambda->bar;
        cout << temp[0] <<" "<< ml->lambda->bar << "\n";
        ss.clear();

        getline( myfile, line );
        ss << line;
        ss >> temp[0] >> temp[1] >> ml->n_atoms;
        cout << temp[0] <<" "<< ml->n_atoms << "\n\n";
        ss.clear();

        ml->atoms = new int[ml->n_atoms];

        getline( myfile, line );
        ss << line;
        ss >> temp[0] >> temp[1];
        for (int k = 0; k < ml->n_atoms; k++)
        {
            ss >> ml->atoms[k];
            ml->atoms[k] -= 1; // to shift index starting from 0
        }
        ss.clear();
        getline( myfile, line );

        // get charges of atoms in lambda group
        ml->chargeA = new real[ml->n_atoms];
        ml->chargeB = new real[ml->n_atoms];
        for (int kk = 0; kk < ml->n_atoms; kk++)
        {
        	getline( myfile, line );
        	ss << line;
        	ss >> temp[0] >> ml->chargeA[kk] >> ml->chargeB[kk];
            cout <<  temp[0] << " " << ml->chargeA[kk] << " " << ml->chargeB[kk] << "\n";
            ss.clear();
        }
        cout << "\n";
		getline( myfile, line ); // empty line

        pos++;
        // open output file for this lambda
        char fname_lambda[255];
        sprintf(fname_lambda, "lambda_%d.dat", (pos));
        ml->out = fopen(fname_lambda, "w");
        out_T = fopen("temperature.dat", "w");

        // initialize lambda structure
        ml->lambda->x    = ml->lambda->x0;
        ml->lambda->T    = cphmd_gen->T_lambda;
        ml->lambda->tau  = cphmd_gen->tau;
        ml->lambda->m    = cphmd_gen->m_lambda;
        ml->lambda->dvdl = 0.0;

        if(cphmd_gen->m_lambda != 0.0)
        {
            // initial velocity of lambda particle
            gmx::UniformRealDistribution<real> uniformDist;
            static int start_seed = 0;
            start_seed = static_cast<int>(gmx::makeRandomSeed());
            gmx::DefaultRandomEngine Rng(start_seed);
            double sigma = sqrt(ml->lambda->T * BOLTZ / ml->lambda->m);
            ml->lambda->v = gaussdist(&Rng, sigma); /* random start velocity to lambda */
            ml->lambda->v_old = 0.0;

            if (ml->residue_name.compare("BUF") == 0)
            {
            	ml->lambda->v = 0.0;
            	ml->lambda->v_old = 0.0;
            }
        }
        else {
        	ml->lambda->v = 0.0;
        	ml->lambda->v_old = 0.0;
        }

        // initialize double well potential
        ml->lambda->lambda_dwp = new real[15];
        init_lambda_dwp(ml->lambda->lambda_dwp, ml->lambda->bar);

        // initial sum of lambdas for charge constraint + buffer mass
		
		for (int j=0; j<cphmd_gen->n_constrained_lambdas; j++ )
		{
		    if (pos==cphmd_gen->constrained_lambdas[j])
		    {
		        if (ml->residue_name.compare("BUF") == 0)
		        {
		        	cphmd_gen->lambda_tot_init = cphmd_gen->lambda_tot_init + cphmd_gen->n_buf*ml->lambda->x0;
		        }
				else
				{
					cphmd_gen->lambda_tot_init = cphmd_gen->lambda_tot_init + ml->lambda->x0;
				}
		    }	
        }

        ml->next = new multi_lambda; // initialize a new multi-lambda structure for the second lambda
        ml       = ml->next;
    }
    
    fprintf(stderr, "lambda_tot_init = %f\n", cphmd_gen->lambda_tot_init);
    fprintf(stderr, "\n");
}

/*  */
void compute_forces(struct cphmd_general *cphmd_gen, struct multi_lambda *ml, real dvdl, real lambdas[], int l_index)
{
    //real  dvdl = 0.0;
    real  corr_ph   = 0.0;
    real  corr_dvdl = 0.0;
    real  corr_dwp  = 0.0;

    real  ref_pka     = cphmd_gen->pKa_values[ml->residue_index]; // reference pKa for current lambda group
    int   n_coeffs    = cphmd_gen->n_coeffs[ml->residue_index];   // number of dvdl coefficicients for current lambda
    real *dvdl_coeffs = new real[n_coeffs];

    // put dvdl coefficients of current lambda group to an array
    for (int k = 0; k < n_coeffs; k++)
    {
        dvdl_coeffs[k] = cphmd_gen->dvdl_coeffs[ml->residue_index*MAX_N_DVDL_COEFFS+k];
    }
	
    // force from electrostatic potential of lambda atoms
    //for (int k = 0; k < ml->n_atoms; k++)
    //{
    //    dvdl += pot[ml->atoms[k]]*(ml->chargeB[k] - ml->chargeA[k]);
    //}
	// dvdl array now has the force so no need to loop over atoms (lambdaAtoms array can be ordered)
	
    // (1) multiple states with three states: two dimensional reference dV/dl
    if(cphmd_gen->multistate_constraint.compare("yes") == 0 && ml->group_number <= cphmd_gen->n_multigroups && cphmd_gen->n_states[ml->group_number-1] > 2)
    {
        int c = 0;
        real l_array[cphmd_gen->n_states[ml->group_number-1]-1];
		
        // save two other lambda coordinates for this lambda dV/dl
        for (int i = 0; i < cphmd_gen->nr_lg; i++)
        {
	        if(i==l_index) continue;
	        l_array[c] = lambdas[i];
	        c += 1;
        }
		
		real *dvdl_coeffs_other = new real[n_coeffs];
		// put other dvdl coefficients of current lambda group to an array
		for (int k = 0; k < n_coeffs; k++)
		{
		    dvdl_coeffs_other[k] = cphmd_gen->dvdl_coeffs_other[ml->residue_index*MAX_N_DVDL_COEFFS+k];
		}
			
		// check which dV/dl to pick, which side of constrained triangle
		if( (lambdas[l_index] + l_array[0]) > (lambdas[l_index] + l_array[1]) ) 
		{
		
		    for (int i = 0; i < n_coeffs; i++)
		    {
		        corr_dvdl += -1.0*dvdl_coeffs[i]*pow(ml->lambda->x, i);
		    }
		}
		else
		{
		    for (int i = 0; i < n_coeffs; i++)
		    {
		        corr_dvdl += -1.0*dvdl_coeffs_other[i]*pow(ml->lambda->x, i);
		    }
		}

    }
	// (2) just two states (multiple states constraint can also be on)
    else
    {   
	    for (int i = 0; i < n_coeffs; i++)
	    {
	        corr_dvdl += -1.0*dvdl_coeffs[i]*pow(ml->lambda->x, i);
	    } 
    }
	
    // effect of external pH
    corr_ph += BOLTZ*cphmd_gen->T_lambda*std::log(10.0)*(ref_pka - cphmd_gen->ph_value);

    // effect of double barrier potential - taken from previous constant ph version 5.1
    real xx = ml->lambda->x;
    real s   = ml->lambda->lambda_dwp[0];
    real f   = ml->lambda->lambda_dwp[1];
    real st  = ml->lambda->lambda_dwp[4];
    real m   = ml->lambda->lambda_dwp[5];
    real a   = ml->lambda->lambda_dwp[6];
    real b   = ml->lambda->lambda_dwp[7];
    real k   = ml->lambda->lambda_dwp[8];
    real d   = ml->lambda->lambda_dwp[9];

    corr_dwp = -k*((-1.0*(xx-1.0-b)*expf(-0.5*(xx-1.0-b)*(xx-1.0-b)/(a*a))/(a*a))
               +(-1.0*(xx+b)*expf(-0.5*(xx+b)*(xx+b)/(a*a))/(a*a)))
               -d*(xx-0.5)*expf(-(xx-0.5)*(xx-0.5)/(2*s*s))/(s*s)
               +f*0.5*(2*st*expf(-(st*st)*(xx-1.0-m)*(xx-1.0-m))/std::sqrt(M_PI)
               -2*st*expf(-(st*st)*(xx+m)*(xx+m))/std::sqrt(M_PI));

    // for dvdl output
    ml->lambda->dvdl_dwp = corr_dwp;
    ml->lambda->dvdl_ref = corr_dvdl;
    ml->lambda->dvdl_ph = corr_ph;
    ml->lambda->dvdl_pot = dvdl;

    // If mass zero then don't update -> reference free energy calculation
    if(cphmd_gen->m_lambda != 0.0)
    {
    	dvdl += corr_ph + corr_dvdl + corr_dwp;
    }
    ml->lambda->dvdl = dvdl;
}

/*
 *  Update lambda and velocity using normal leap-frog
 *  - needs thermostat (collective v-rescale)
 */
void updateLambda(const t_inputrec *ir, struct multi_lambda *ml, struct cphmd_general *cphmd_gen)
{
	if(cphmd_gen->m_lambda != 0.0)
	{
		real inverseMass = 1.0/ml->lambda->m;

		ml->lambda->v = ml->lambda->v + inverseMass*(-1.0)*ml->lambda->dvdl*ir->delta_t;
        //fprintf(stderr, "v updateLambda = %g\n", ml->lambda->v);
		ml->lambda->ekin = 0.5*ml->lambda->m*ml->lambda->v*ml->lambda->v;
		ml->lambda->T = 2.0*(ml->lambda->ekin)/BOLTZ;
		ml->lambda->x_old = ml->lambda->x;
		ml->lambda->x = ml->lambda->x_old + ml->lambda->v*ir->delta_t;
	}
}

/*
 *  Update lambda and velocity using stochastic (Langevin) dynamics
 *  - modified from gromacs langevin function doSDUpdateGeneral() in update.cpp
 *  - phases 1,2 take into account the charge constraint
 */
void updateLambdaLD(const t_inputrec *ir, struct multi_lambda *ml, struct cphmd_general *cphmd_gen, real T_ref, real gamma, int phase)
{

    if(cphmd_gen->m_lambda != 0.0)
    {
    	real inverseMass = 1.0/ml->lambda->m;

    	// compare to updateType == SDUpdate::ForcesOnly in update.cpp
        if( phase == 1 )
        {
        	ml->lambda->x_old_old = ml->lambda->x;
        	ml->lambda->v_old_old = ml->lambda->v;

    	    ml->lambda->v = ml->lambda->v + inverseMass*(-1.0)*ml->lambda->dvdl*ir->delta_t;
    	    ml->lambda->x = ml->lambda->x + ml->lambda->v*ir->delta_t;

    	    ml->lambda->x_old = ml->lambda->x;
    	    ml->lambda->v_old = ml->lambda->v;
        }

        // compare to updateType == SDUpdate::FrictionAndNoiseOnly in update.cpp
        if( phase == 2 )
        {
	        int seed = static_cast<int>(gmx::makeRandomSeed());
	        // copied from update.cpp
	        gmx::ThreeFry2x64<0> rng(seed, gmx::RandomDomain::UpdateCoordinates);
	        gmx::TabulatedNormalDistribution<real, 14> dist;

	        real deltav, f;
	        real kT = BOLTZ*T_ref;

	        real inverseSqrMass = std::sqrt(inverseMass);

	        f = 1.0-std::exp(-1.0*gamma*ir->delta_t);

	        ml->lambda->v_old_old = ml->lambda->v;

	        deltav = (1.0-f)*ml->lambda->v + inverseSqrMass*std::sqrt(kT*f*(2.0-f))*dist(rng);
	        ml->lambda->x = ml->lambda->x + 0.5*(deltav - ml->lambda->v)*ir->delta_t;

	        ml->lambda->x_old_old = ml->lambda->x;

	        //ml->lambda->v_old = 0.5*(deltav - ml->lambda->v);
	        ml->lambda->v_old = deltav;
	        ml->lambda->x_old = ml->lambda->x;
        }

        // in update.cpp the "normal" update without constraint, add friction and noise
        // here called if charge constraint is not on
        if( phase == 3 )
        {
            int seed = static_cast<int>(gmx::makeRandomSeed());
        	// copied from update.cpp
        	gmx::ThreeFry2x64<0> rng(seed, gmx::RandomDomain::UpdateCoordinates);
        	gmx::TabulatedNormalDistribution<real, 14> dist;

        	real deltav, f;
        	real kT = BOLTZ*T_ref;

        	real inverseSqrMass = std::sqrt(inverseMass);

        	f = 1.0-std::exp(-1.0*gamma*ir->delta_t);

        	ml->lambda->v = ml->lambda->v + inverseMass*(-1.0)*ml->lambda->dvdl*ir->delta_t;
        	deltav = -1.0*f*ml->lambda->v + inverseSqrMass*std::sqrt(kT*f*(2.0-f))*dist(rng);

        	ml->lambda->x = ml->lambda->x + (ml->lambda->v + 0.5*deltav)*ir->delta_t;
        	ml->lambda->v = ml->lambda->v + deltav;

            ml->lambda->ekin = 0.5*ml->lambda->m*ml->lambda->v*ml->lambda->v;
            ml->lambda->T = 2.0*(ml->lambda->ekin)/BOLTZ;
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
void do_constraints(struct cphmd_general *cphmd_gen, struct multi_lambda *ml, real dt)
{	
	multi_lambda *tmp = new multi_lambda;
	tmp = ml;
	
	// (1) multiple states current lambda
	
	real l_multi_init = 1.0;
	real l_multi[cphmd_gen->n_multigroups];
	real deviation [cphmd_gen->n_multigroups];
	real max_deviation;
	
	if (cphmd_gen->multistate_constraint.compare("no")==0) 
	{
		max_deviation = 0.0;
	}
	else
	{ 
	    // compute current sum of lambdas for all multistates groups
	    for (int group = 1; group <= cphmd_gen->n_multigroups; group++)
	    {
	        for (int i = 0; i < cphmd_gen->nr_lg; i++)
	        {
		        // check if lambda groups belongs to certain multiple states group
		        if (ml->group_number == group)
    		    {
	    		    l_multi[group-1] = l_multi[group-1] + ml->lambda->x;
		        }
		        ml = ml->next;
		    }
			ml = tmp;
	    }

	    for(int i=0; i<cphmd_gen->n_multigroups; i++) 
	    {
		    deviation[i] = abs(l_multi[i]-l_multi_init);
	    }
	    max_deviation = *std::max_element(deviation, deviation + cphmd_gen->n_multigroups-1); // largest deviation 
    }

	// (2) charge constraint current total lambda
	
	real l_cc_init = cphmd_gen->lambda_tot_init;
	real l_cc = 0;
	
    if (cphmd_gen->charge_constraint.compare("no")==0) 
    {
	    l_cc = l_cc_init;
    }
    else
    {
	    // compute current sum of all lambdas and collective buffer
    	for (int i = 0; i<cphmd_gen->nr_lg; i++)
	    {
		    for (int j=0; j<cphmd_gen->n_constrained_lambdas; j++ )
		    {
		        if (i==cphmd_gen->constrained_lambdas[j]-1)
		        {
		            if (ml->residue_name.compare("BUF") == 0)
		            {
		            	l_cc = l_cc + cphmd_gen->n_buf*ml->lambda->x;
		            }
			    	else
				    {
					    l_cc = l_cc + ml->lambda->x;
				    }
		        }	
	    	}
		    ml = ml->next;
    	}
    	ml = tmp;
    }	
	
	// iterate until all multiple states and charge constraint are close to desired
	int it = 0;
	real multiplier;
	real sum_dt2_per_mass;
	int counter;  // how many lambda groups in multistate group
	
	while ( (it<10000 && max_deviation > 0.00001) || (it<10000 && abs(l_cc - l_cc_init) > 0.00001) ) 
	{
		// multiple states constraint for all multistate groups
		
		if(cphmd_gen->multistate_constraint.compare("yes")==0)
		{
			// loop over all multistate groups
    	    for (int group = 1; group <= cphmd_gen->n_multigroups; group++)
    	    {
				multiplier = 0;  // Lagrange multiplier
				l_multi[group-1] = 0;     // sum of multigroup lambdas
				
				// sum current multistate group lambdas
				counter = 0;
				for (int i = 0; i < cphmd_gen->nr_lg; i++)
				{
				    // check if lambda groups belongs to this multiple states group
				    if (ml->group_number == group)
				    {
				        l_multi[group-1] = l_multi[group-1] + ml->lambda->x;
						counter = counter + 1;
				    }
				    ml = ml->next;
				}
				ml = tmp;
				
				sum_dt2_per_mass = counter*dt*dt/cphmd_gen->m_lambda;
				multiplier = (l_multi[group-1] - l_multi_init)/sum_dt2_per_mass;
				l_multi[group-1] = 0;
				
				// linear constraint for one group converges in one step
				for (int i = 0; i < cphmd_gen->nr_lg; i++)
				{
					if (ml->group_number == group)
					{
					    ml->lambda->x =  ml->lambda->x - multiplier*dt*dt/ml->lambda->m;
						l_multi[group-1] = l_multi[group-1] + ml->lambda->x;
					} 
					ml = ml->next;
				}
				ml = tmp;  
    	    }
	
			for(int i=0; i<cphmd_gen->n_multigroups; i++) 
			{
				deviation[i] = abs(l_multi[i]-l_multi_init);
				//fprintf(stderr,"group %i deviation %g \n",i+1,deviation[i]);
			}
			max_deviation = *std::max_element(deviation, deviation + cphmd_gen->n_multigroups-1); // largest deviation from l_multi_init = 1
			//fprintf(stderr,"max deviation %g \n\n",max_deviation);
		}
		
		// charge constraint for total sum of lambdas
		
		if(cphmd_gen->charge_constraint.compare("yes")==0)
		{
			sum_dt2_per_mass = (cphmd_gen->n_constrained_lambdas-1)*dt*dt/cphmd_gen->m_lambda + cphmd_gen->n_buf*cphmd_gen->n_buf*dt*dt/cphmd_gen->m_buf;
			multiplier = 0;  // Lagrange multiplier
			l_cc = 0;
			
			// sum all lambdas (and collective buffer)
			for (int i = 1; i<=cphmd_gen->nr_lg; i++)
			{
				for (int j=0; j<cphmd_gen->n_constrained_lambdas; j++ )
				{
				    if (i==cphmd_gen->constrained_lambdas[j])
				    {
				        if (ml->residue_name.compare("BUF") == 0)
				        {	
				        	l_cc = l_cc + cphmd_gen->n_buf*ml->lambda->x;
				        }
						else
						{
							l_cc = l_cc + ml->lambda->x;
						}
				    }	
				}
				ml = ml->next;
			}
			ml = tmp;

			multiplier = (l_cc - l_cc_init)/sum_dt2_per_mass;
			l_cc = 0;
			
			// linear constraint converges in one step
			for (int i = 1; i<=cphmd_gen->nr_lg; i++)
			{
				for (int j=0; j<cphmd_gen->n_constrained_lambdas; j++ )
				{
				    if (i==cphmd_gen->constrained_lambdas[j])
				    {
				        if (ml->residue_name.compare("BUF") == 0)
				        {
					        ml->lambda->x =  ml->lambda->x - cphmd_gen->n_buf*multiplier*dt*dt/ml->lambda->m;
					        l_cc = l_cc + cphmd_gen->n_buf*ml->lambda->x;
				        }
				        else
				        {
					        ml->lambda->x =  ml->lambda->x - multiplier*dt*dt/ml->lambda->m;
				            l_cc = l_cc + ml->lambda->x;
				        }
			        }
		        }
				ml = ml->next;
			 }
		     ml = tmp;
		}	
		it += 1;
	}
	
	// advance lambda velocities, kinetic energies, temperatures
	// (also for non-constrained groups these now computed again but should not matter)
	
	for (int i = 0; i < cphmd_gen->nr_lg; i++)
	{
		ml->lambda->v = (ml->lambda->x - ml->lambda->x_old)/dt;
		//fprintf(stderr,"x = %g x_old = %g \n",ml->lambda->x, ml->lambda->x_old);
		//fprintf(stderr,"v = %g \n",ml->lambda->v);
		ml->lambda->ekin = 0.5*ml->lambda->m*ml->lambda->v*ml->lambda->v;
		ml->lambda->T = 2.0*(ml->lambda->ekin)/BOLTZ;
		ml = ml->next;
	}
	ml = tmp; 
}

/*
 * Constraint for multiple sites 
 */
void do_multiple_states_constraint(struct cphmd_general *cphmd_gen, struct multi_lambda *ml, real dt, int group)
{
	real multiplier = 0;  // Lagrange multiplier
	real lambda_tot = 0;  // total sum of lambdas
	int counter = 0;  // how many lambda groups in current multistate group
	multi_lambda *tmp = new multi_lambda;
	tmp = ml;
	
	// sum current lambda values
	for (int i = 0; i < cphmd_gen->nr_lg; i++)
	{
		// check if lambda groups belongs to certain multiple states group
		if (ml->group_number == group)
		{
			lambda_tot = lambda_tot + ml->lambda->x;
			counter = counter + 1;
		}
		ml = ml->next;
	}
	ml = tmp;
	
	real sum_dt2_per_mass = counter*dt*dt/cphmd_gen->m_lambda;

	// linear constraint converges in one step
	real lambda_tot_goal = 1.0;
	multiplier = (lambda_tot - lambda_tot_goal)/sum_dt2_per_mass;
	for (int i = 0; i < cphmd_gen->nr_lg; i++)
	{
		if (ml->group_number == group)
		{
		    ml->lambda->x =  ml->lambda->x - multiplier*dt*dt/ml->lambda->m;
		}
		ml = ml->next;
	}
    ml = tmp;

    for (int i = 0; i < cphmd_gen->nr_lg; i++)
	{
		if (ml->group_number == group)
		{
			ml->lambda->v = (ml->lambda->x - ml->lambda->x_old)/dt;
			ml->lambda->ekin = 0.5*ml->lambda->m*ml->lambda->v*ml->lambda->v;
			ml->lambda->T = 2.0*(ml->lambda->ekin)/BOLTZ;
		}
		ml = ml->next;
	}
	ml = tmp;

}

void do_charge_constraint(struct cphmd_general *cphmd_gen, struct multi_lambda *ml, real sum_dt2_per_mass, real dt, int phase)
{
	real multiplier = 0;  // Lagrange multiplier
	real lambda_tot = 0;  // total sum of lambdas
	real lambda_tot_init = cphmd_gen->lambda_tot_init; //;
	multi_lambda *tmp = new multi_lambda;
	tmp = ml;

	// sum current lambda values
	for (int i = 0; i < cphmd_gen->nr_lg; i++)
	{
		if (i == (cphmd_gen->nr_lg - 1))
		{
		    lambda_tot = lambda_tot + cphmd_gen->n_buf*ml->lambda->x;
		}
		else
		{
			lambda_tot = lambda_tot + ml->lambda->x;
		}
		ml = ml->next;
	}
	ml = tmp;

	// linear constraint converges in one step
	multiplier = (lambda_tot - lambda_tot_init)/sum_dt2_per_mass;
	for (int i = 0; i < cphmd_gen->nr_lg; i++)
	{
		if (i == (cphmd_gen->nr_lg - 1))
		{
			ml->lambda->x =  ml->lambda->x - cphmd_gen->n_buf*multiplier*dt*dt/ml->lambda->m;
		}
		else
		{
			ml->lambda->x =  ml->lambda->x - multiplier*dt*dt/ml->lambda->m;
		}
		ml = ml->next;
	    }
		ml = tmp;

	// phase 0 for v-rescale thermostat
	if (phase==0)
		{
			for (int i = 0; i < cphmd_gen->nr_lg; i++)
			{
				ml->lambda->v = (ml->lambda->x - ml->lambda->x_old)/dt;
                //fprintf(stderr, "v constraint = %g\n", ml->lambda->v);
				ml->lambda->ekin = 0.5*ml->lambda->m*ml->lambda->v*ml->lambda->v;
				ml->lambda->T = 2.0*(ml->lambda->ekin)/BOLTZ;
			    ml = ml->next;
			}
			ml = tmp;
		}

	// phases 1,2 for langevin charge constraint
	if (phase==1)
	{
		for (int i = 0; i < cphmd_gen->nr_lg; i++)
		{
			ml->lambda->v = (ml->lambda->x - ml->lambda->x_old_old)/dt;
		    ml = ml->next;
		}
		ml = tmp;
	}
	if (phase==2)
	{
		for (int i = 0; i < cphmd_gen->nr_lg; i++)
		{
			ml->lambda->v = ml->lambda->v_old_old + 2.0*(ml->lambda->x - ml->lambda->x_old_old)/dt;

			ml->lambda->T = 2.0*(ml->lambda->ekin)/BOLTZ;
		    ml = ml->next;
		}
	    ml = tmp;
	}
}

/*  Perform T coupling with v-rescale thermostat
 *  - couple all lambdas together
 *  - Nf = number of degrees of freedom = N - 1
 *  - according to Bussi et al JCP (2007) appendix
 *
 *  Returns the change in kinetic energy
 */
real tcouple_vrescale_collective(struct cphmd_general *cphmd_gen, struct multi_lambda *ml, real Tref, real tau, real dt, int64_t step, int64_t seed)
{
    real alpha, factor, sum_r2, r, r1, Ekin_total,Ekin_ref;
    int Nf;
    gmx::ThreeFry2x64<64> rng(seed, gmx::RandomDomain::Thermostat);

    gmx::NormalDistribution<real> normalDist;
    multi_lambda *tmp = new multi_lambda;
    tmp = ml;

    factor = exp(-1.0*dt/tau);

    // number of degrees of freedom (less if constraint on)
    Nf = cphmd_gen->nr_lg;
    if(cphmd_gen->charge_constraint.compare("yes")==0)
    {
    	Nf = Nf - 1;
    }
    if(cphmd_gen->multistate_constraint.compare("yes")==0)
    {
    	Nf = Nf - cphmd_gen->n_multigroups;
    }

    // total kinetic energy
    Ekin_total = 0;
    for (int i = 0; i < cphmd_gen->nr_lg; i++)
    {
        Ekin_total = Ekin_total + 0.5*ml->lambda->m*ml->lambda->v*ml->lambda->v;
    	ml = ml->next;
    }
    ml = tmp;

    // reference kinetic energy from desired temperature
    Ekin_ref = Nf*0.5*Tref*BOLTZ;

    // first random number
    rng.restart(step, 0);
    r1 = normalDist(rng);

    // sum of gaussian squared random numbers
    sum_r2 = 0;
    for(int j = 1; j < Nf; j++)
    {
    	r = normalDist(rng);
    	sum_r2 = sum_r2 + r*r;
    }
    real Ekin_new;
    Ekin_new = Ekin_total + (1.0-factor)*(Ekin_ref*(r1*r1+sum_r2)/Nf - Ekin_total) + 2.0*r1*std::sqrt(Ekin_ref*Ekin_total/Nf*(1.0-factor)*factor);

    // Analytically Ek_new>=0, but we check for rounding errors (from gromacs coupling.cpp)
    if (Ekin_new <= 0)
    {
        alpha = 0.0;
    }
    else
    {
    	alpha = std::sqrt(Ekin_new/Ekin_total);
    }

    // scale kinetic energies and velocities with alpha
    for (int i = 0; i < cphmd_gen->nr_lg; i++)
    {
		ml->lambda->v = alpha*ml->lambda->v;
        //fprintf(stderr, "v v-rescale = %g\n", ml->lambda->v);
		// new positions
		ml->lambda->x = ml->lambda->x_old + ml->lambda->v*dt;
		ml->lambda->ekin = alpha*alpha*ml->lambda->ekin;
		ml->lambda->T = 2.0*(ml->lambda->ekin)/BOLTZ;

    	ml = ml->next;
    }
    ml = tmp;

    return Ekin_new - Ekin_total;
}

/* Initialize double well barrier potential  
   Previous work taken from constant pH version 5.1 
*/
void init_lambda_dwp(real *lambda_dwp, real barrier)
{

    /* fields of array lambda_dwp are filled 08.08.2014 */
    FILE* fw1;
    real  s = 0.3;
    real  e1, e10, st, m;
    real  sig0, sig, x, x0;
    real  eps = 0.005;
    real  dg;
    real  v, vmin, vmin0, dv, iii;
    real  totx, totp, totxp;
    real  flagfor;
    int   ready;
    real  a, b, k, d;
    int   iter;
    int   max_iter;
    b        = -0.1;
    a        = 0.05;
    iter     = 0;
    max_iter = 10000;

    fw1 = fopen("lambda_dwp.dat", "w");
    iii = barrier - 1.0;
    sig0  = 0.02;
    sig   = sig0;
    x0    = 0;
    vmin0 = 10.0;

    /* Note that the height of outer walls will be f+barrier/2 */

    // for f=50 follows
    //  real f=50.0;
    //  e1=1.45222;         /*e1=InverseErf[1-2.0/f]*/
    // e10=0.595116;        /*e10=InverseErf[1-20.0/f]*/

    //real f=100.0;
    //e1=1.64498;           /*e1=InverseErf[1-2.0/f]*/
    //e10=0.906194;         /*e10=InverseErf[1-20.0/f]*/

    // For f=1000 follows
    real f = 1000.0;
    e1  = 2.185;
    e10 = 1.645;

    lambda_dwp[0] = s;
    lambda_dwp[1] = f;
    lambda_dwp[2] = e1;
    lambda_dwp[3] = e10;

    st = (e1 - e10) / (2.0 * sig0);
    m  = 2.0 * sig0 * (2.0 * e1 - e10) / (e1 - e10);

    lambda_dwp[4] = st;
    lambda_dwp[5] = m;

    flagfor = iii;

    for (iii = flagfor; iii <= flagfor; iii = iii + 1)
    {
        dg   = 1.0 + iii;
        sig0 = 0.02;
        sig  = sig0;
        b    = -0.1;
        a    = 0.05;
        k    = dg / 2.0;
        d    = dg / 2.0;

        /*  v=-k*(expf(-((x-1.0-b)*(x-1.0-b)/(2*a*a)))+expf(-(x+b)*(x+b)/(2*a*a)))+
         *       d*expf(-(x-0.5)*(x-0.5)/(2*s*s))+
         *            f*0.5*((1-std::erf(st*(x+m)))+(1.0+std::erf(st*(x-1.0-m))));
         *                this is the barrier function */

        /* correct for the fact that the two minima are shallower than k=dg/2 due to the tails of the
         *    central gaussian and the erf functions
         *       (basically bring the minima to the same value of the two side gaussians
         *          before the central gaussian was added) */

        vmin = vmin0;
        for (x = -0.1; x <= 0.2; x = x + 0.001)
        {
            v = -k * (std::exp(-((x - 1.0 - b) * (x - 1.0 - b) / (2 * a * a)))
                      + std::exp(-(x + b) * (x + b) / (2 * a * a)))
                + d * std::exp(-(x - 0.5) * (x - 0.5) / (2 * s * s))
                + f * 0.5
                * ((1 - std::erf(st * (x + m)))
                   + (1.0 + std::erf(st * (x - 1.0 - m))));
            if (v < vmin)
            {
                vmin = v;
            }
        }
        k = k + dg / 2.0 + vmin;

        /* adjust location minima and width boltzmann distribution therein to the target values (0,1 and sig0) */

        ready = 0;
        while (ready == 0)   /*while which correspond to repeat begin*/

        {
            vmin = vmin0;
            for (x = -0.1; x < 0.2; x = x + 0.001)
            {
                v = -k * (std::exp(-((x - 1.0 - b) * (x - 1.0 - b) / (2 * a * a)))
                          + std::exp(-(x + b) * (x + b) / (2 * a * a)))
                    + d * std::exp(-(x - 0.5) * (x - 0.5) / (2 * s * s))
                    + f * 0.5
                    * ((1 - std::erf(st * (x + m)))
                       + (1.0 + std::erf(st * (x - 1.0 - m))));
                if (v < vmin)
                {
                    vmin = v;
                }
            }
            k = k + dg / 2.0 + vmin;
            /*v=-k*(std::exp(-((x-1.0-b)*(x-1.0-b)/(2*a*a)))+std::exp(-(x+b)*(x+b)/(2*a*a)))+
             *         d*std::exp(-(x-0.5)*(x-0.5)/(2*s*s))+
             *                 f*0.5*((1-std::erf(st*(x+m)))+(1.0+std::erf(st*(x-1.0-m))));*/

            x    = 0.5;
            totp = 0;
            totx = 0;
            while (x > -0.2)
            {
                v = -k * (std::exp(-((x - 1.0 - b) * (x - 1.0 - b) / (2 * a * a)))
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
                x = x - 0.001;
            }
            x0 = totx / totp;

            x     = 0.5;
            totxp = 0;
            while (x > -0.2)
            {
                v = -k * (std::exp(-((x - 1.0 - b) * (x - 1.0 - b) / (2 * a * a)))
                          + std::exp(-(x + b) * (x + b) / (2 * a * a)))
                    + d * std::exp(-(x - 0.5) * (x - 0.5) / (2 * s * s))
                    + f * 0.5
                    * ((1 - std::erf(st * (x + m)))
                       + (1.0 + std::erf(st * (x - 1.0 - m))));
                if (v <= 0)
                {
                    totxp = totxp + ((x - x0) * (x - x0) * std::exp(-v));
                }
                x = x - 0.001;
            }
            sig = std::sqrt(totxp / totp);
            b   = b + 0.01 * x0;
            a   = a / (1.0 + 0.01 * (sig - sig0) / sig0);

            for (x = -0.2; x <= 1.2; x = x + 0.001)
            {
                v = -k * (std::exp(-((x - 1.0 - b) * (x - 1.0 - b) / (2 * a * a)))
                          + std::exp(-(x + b) * (x + b) / (2 * a * a)))
                    + d * std::exp(-(x - 0.5) * (x - 0.5) / (2 * s * s))
                    + f * 0.5
                    * ((1 - std::erf(st * (x + m)))
                       + (1.0 + std::erf(st * (x - 1.0 - m))));
            }

            if ((std::abs(x0) <= eps) && (std::abs(sig - sig0) / sig0 <= eps))
            {
                ready = 1;
                for (x = -1.0; x <= 2.0; x = x + 0.001)
                {
                    v = -k * (std::exp(-((x - 1.0 - b) * (x - 1.0 - b)
                                         / (2 * a * a)))
                              + std::exp(-(x + b) * (x + b) / (2 * a * a)))
                        + d * std::exp(-(x - 0.5) * (x - 0.5) / (2 * s * s))
                        + f * 0.5
                        * ((1 - std::erf(st * (x + m)))
                           + (1.0 + std::erf(st * (x - 1.0 - m))));

                    dv =
                        -k
                        * ((-1.0 * (x - 1.0 - b)
                            * std::exp(
                                    -0.5 * (x - 1.0 - b)
                                    * (x - 1.0 - b)
                                    / (a * a)) / (a * a))
                           + (-1.0 * (x + b)
                              * std::exp(
                                      -0.5 * (x + b)
                                      * (x + b)
                                      / (a * a))
                              / (a * a)))
                        - d * (x - 0.5)
                        * std::exp(
                                -(x - 0.5) * (x - 0.5)
                                / (2 * s * s))
                        / (s * s)
                        + f * 0.5
                        * (2 * st
                           * std::exp(
                                   -(st * st)
                                   * (x - 1.0
                                      - m)
                                   * (x - 1.0
                                      - m))
                           / std::sqrt(M_PI)
                           - 2 * st
                           * std::exp(
                                   -(st * st)
                                   * (x
                                      + m)
                                   * (x
                                      + m))
                           / std::sqrt(M_PI));
                    fprintf(fw1, "%f %f\n", x, v);
                }
                lambda_dwp[6] = a;
                lambda_dwp[7] = b;
                lambda_dwp[8] = k;
                lambda_dwp[9] = d;
                /* fprintf(stderr,
                        "Parameters double well potential (applied only if double well potential is on):\n");
                fprintf(stderr, "a = %f b = %f k = %f d = %f dg = %f \n", a, b,
                        k, d, dg);
                fprintf(stderr, "m = %f st = %f \n", m, st);
                */
            }

            if (iii < -0.55)
            {
                k     = 0;
                d     = 0;
                ready = 1;
                for (x = -1.0; x <= 2.0; x = x + 0.001)
                {
                    v = -k * (std::exp(-((x - 1.0 - b) * (x - 1.0 - b)
                                         / (2 * a * a)))
                              + std::exp(-(x + b) * (x + b) / (2 * a * a)))
                        + d * std::exp(-(x - 0.5) * (x - 0.5) / (2 * s * s))
                        + f * 0.5
                        * ((1 - std::erf(st * (x + m)))
                           + (1.0 + std::erf(st * (x - 1.0 - m))));

                    dv = -k * ((-1.0 * (x - 1.0 - b)
                                * std::exp(
                                        -0.5 * (x - 1.0 - b)
                                        * (x - 1.0 - b)
                                        / (a * a)) / (a * a))
                               + (-1.0 * (x + b)
                                  * std::exp(
                                          -0.5 * (x + b)
                                          * (x + b)
                                          / (a * a))
                                  / (a * a)))
                        - d * (x - 0.5)
                        * std::exp(
                                -(x - 0.5) * (x - 0.5)
                                / (2 * s * s))
                        / (s * s)
                        + f * 0.5
                        * (2 * st
                           * std::exp(
                                   -(st * st)
                                   * (x - 1.0
                                      - m)
                                   * (x - 1.0
                                      - m))
                           / std::sqrt(M_PI)
                           - 2 * st
                           * std::exp(
                                   -(st * st)
                                   * (x
                                      + m)
                                   * (x
                                      + m))
                           / std::sqrt(M_PI));
                    fprintf(fw1, "%f %f\n", x, v);
                }
                lambda_dwp[6] = a;
                lambda_dwp[7] = b;
                lambda_dwp[8] = k;
                lambda_dwp[9] = d;
                if (barrier != 0)
                {
                    fprintf(stderr,
                            "Warning: Barrier changed from %f to zero\n",
                            barrier);
                }
                fprintf(stderr, "Parameters double well potential:\n");
                fprintf(stderr, "a = %f b = %f k = %f d = %f dg = %f \n", a, b,
                        k, d, dg);
                fprintf(stderr, "m = %f st = %f \n", m, st);
				fprintf(stderr, "\n");
            }

            iter++;

            if (iter > max_iter)
            {
                ready = 1;
                for (x = -1.0; x <= 2.0; x = x + 0.001)
                {
                    v = -k * (std::exp(
                                      -((x - 1.0 - b) * (x - 1.0 - b)
                                        / (2 * a * a)))
                              + std::exp(-(x + b) * (x + b) / (2 * a * a)))
                        + d * std::exp(-(x - 0.5) * (x - 0.5) / (2 * s * s))
                        + f * 0.5
                        * ((1 - std::erf(st * (x + m)))
                           + (1.0 + std::erf(st * (x - 1.0 - m))));
                    
                    dv = -k * ((-1.0 * (x - 1.0 - b)
                                * std::exp(
                                        -0.5 * (x - 1.0 - b)
                                        * (x - 1.0 - b)
                                        / (a * a)) / (a * a))
                               + (-1.0 * (x + b)
                                  * std::exp(
                                          -0.5 * (x + b)
                                          * (x + b)
                                          / (a * a))
                                  / (a * a)))
                        - d * (x - 0.5)
                        * std::exp(
                                -(x - 0.5) * (x - 0.5)
                                / (2 * s * s))
                        / (s * s)
                        + f * 0.5
                        * (2 * st
                           * std::exp(
                                   -(st * st)
                                   * (x - 1.0
                                      - m)
                                   * (x - 1.0
                                      - m))
                           / std::sqrt(M_PI)
                           - 2 * st
                           * std::exp(
                                   -(st * st)
                                   * (x
                                      + m)
                                   * (x
                                      + m))
                           / std::sqrt(M_PI));
                    fprintf(fw1, "%f %f %f\n", x, v,dv);
                }
                lambda_dwp[6] = a;
                lambda_dwp[7] = b;
                lambda_dwp[8] = k;
                lambda_dwp[9] = d;
                fprintf(stderr,
                        "Warning: double well potential did not converge to tolerance limit eps=%f in max iter=%d\n",
                        eps, max_iter);
                fprintf(stderr,
                        "Warning: check the shape of the double well potential in lambda_dwp.dat\n");
                fprintf(stderr, "Parameters double well potential:\n");
                fprintf(stderr, "a = %f b = %f k = %f d = %f dg = %f \n", a, b,
                        k, d, dg);
                fprintf(stderr, "m = %f st = %f \n", m, st);
				fprintf(stderr, "\n");
            }

        } /* endwhile of repeat begin */

    }     /* end for */

    fclose(fw1);

}
