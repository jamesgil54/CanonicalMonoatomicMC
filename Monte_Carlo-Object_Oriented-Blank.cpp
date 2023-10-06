/* 
This code runs Monte Carlo simulations for Lennard Jones particles of a single type, in the canonical (NVT) ensemble.

The simulation cell is assumed to be cubic.

R. K. Lindsey (2023)
*/


#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>

#include "MersenneTwiser.h"

using namespace std;

const double  nav     = 6.02*pow(10.0,23);        // Avagadro's number
const double  cm2m    = 1.0/100.0;                // x cm * cm2m = m
const double  A2m     = 1.0/pow(10.0,10);         // x Angstrom * A to m = m    
const double  kB      = 1.380649*pow(10,-23.0);   // Units: J⋅K^−1 

struct xyz
{
    double x;
    double y;
    double z;
};

double get_dist(const xyz & a1, const xyz & a2, const xyz & boxdim, xyz & rij_vec)
{
    /* Write code to calculate the distance vector and scalar distance between two atoms, using the minimum image convention.
    The distance vector is rij_vec - it is passed by reference so it can be modified directly.
    The convention is that the vector should point from atom 1 to atom 2
    
    The function should return the scalar distance.
    */

    double dx = (a2.x-a1.x) - nint((a2.x-a1.x)/boxdim.x) * boxdim.x;
    double dy = (a2.y-a1.y) - nint((a2.y-a1.y)/boxdim.y) * boxdim.y;
    double dz = (a2.z-a1.z) - nint((a2.z-a1.z)/boxdim.z) * boxdim.z;

    rij_vec.x = dx;
    rij_vec.y = dy;
    rij_vec.z = dz;

    double distance = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dy, 2));

    return distance;
}

double update_max_displacement(double fraction_accepted, double boxdim, double max_displacement)
{
    /* Write code to update the maximum displacement based on current acceptance crtieria. 
    The function should return the the new maximum displacement.
    */
    if(fraction_accepted >= 0.5){
        max_displacement =  max_displacement*1.2;
    }

    if(fraction_accepted < 0.5){
        max_displacement = max_displacement*0.8;
    }

    return max_displacement;
}

class system_coordinates
{
    
public: 
    // General system coordinate definitions
    
    int         natoms;     // Number of atoms in the system
    string      atmtyp;     // Atom type (atomic symbol)
    double      density;    // System density in g/cm^3
    double      numden;     // System number density (atoms/Ang^3)
    double      molmass;    // Molar mass of an atom of atmtyp
    xyz         boxdim;     // Simulation box x, y, and z lengths
    vector<xyz> coords;     // Coordinates for all atoms in our system
    
    // Variables used for file i/o
    
    ofstream    trajstream; // Trajectory file - uses the LAMMPS lammpstrj format

    void generate_coords();
    void write_frame(int mcstep);
    
    // Constructor and deconstructor
    
    system_coordinates(int in_natoms, string in_atmtyp, double in_density, double in_molmass);
    ~system_coordinates();
};

system_coordinates::system_coordinates(int in_natoms, string in_atmtyp, double in_density, double in_molmass)
{
    natoms  = in_natoms;
    atmtyp  = in_atmtyp;
    density = in_density;
    molmass = in_molmass;
    
    trajstream.open("MC_traj.lammpstrj");
}
system_coordinates::~system_coordinates()
{
    trajstream.close();
}

void system_coordinates::generate_coords()
 ///////////////////////////////////////////////////////////////////////////////////////////////
 // Generate the system coordinates
 // Generates an initial configuation of atom on a 3D lattice  with a user-specified number of 
 // atoms, at the user-specified density. 
 ///////////////////////////////////////////////////////////////////////////////////////////////
{
    // Determine the box length that will yield the user-specified density. Start by computing the target number density.

    
    numden = density / molmass * nav / pow(cm2m,3.0) * pow(A2m,3.0); // Density in atoms per Ang^3
    
    cout << "# Num. den. (atoms/Ang): " << numden << endl;
    
    boxdim.x = pow(natoms/numden, 1.0/3.0); 
    boxdim.y = boxdim.x;
    boxdim.z = boxdim.x;
    
    cout << "# Box length (x):        " << boxdim.x << endl;
    cout << "# Box length (y):        " << boxdim.y << endl;
    cout << "# Box length (z):        " << boxdim.z << endl;

    // Compute the number of gridpoints to use in each direction (ceiling of the cube root of number of atoms).
    // Set the spacing in each dimension based on the box length and the number of gridpoints+1 (to prevent overlap
    // across the periodic boundary. 
    
    int     ngridpoints = ceil(pow(natoms,1.0/3.0));
    double  spacing     = boxdim.x / (ngridpoints+1);
    
    cout << "# Init. spacing (Ang):   " << spacing << endl;

    xyz     tmp_atom;
    int     added_atoms = 0;
    
    for (int x=0; x<ngridpoints; x++)
    {
        for(int y=0; y<ngridpoints; y++)
        {
            for(int z=0; z<ngridpoints; z++)  
            {
                if(added_atoms >= natoms)
                    continue;
                
                tmp_atom.x = x * spacing;
                tmp_atom.y = y * spacing;
                tmp_atom.z = z * spacing;

                coords.push_back(tmp_atom);
                
                added_atoms++;
            }
        }
    }
                
    
    // Print this initial configuration to a file
    
    write_frame(-1);   
}

void system_coordinates::write_frame(int mcstep)
{
    trajstream << "ITEM: TIMESTEP" << endl;
    trajstream << mcstep << endl;
    trajstream << "ITEM: NUMBER OF ATOMS" << endl;
    trajstream << natoms << endl;
    trajstream << "ITEM: BOX BOUNDS pp pp pp" << endl;
    trajstream << "0 " << boxdim.x << endl;
    trajstream << "0 " << boxdim.y << endl;
    trajstream << "0 " << boxdim.z << endl;
    trajstream << "ITEM: ATOMS id type element xu yu zu" << endl;
    
    for (int i=0; i<natoms; i++)
        trajstream << i << " 1 " << atmtyp << " " << coords[i].x << " "<< coords[i].y << " "  << coords[i].z << endl;
}

class LJ_model
{
public:
    
    ///units for epsilon: e/Kb?
    
    double  sigma;
    double  epsilon;
    double  rcut;
    
    double  term1;   // Placeholder for (sigma/rij)^12
    double  term2;   // Placeholder for (sigma/rij)^6
    
    double  get_eij(double rij);
    void    get_fij(double rij, const xyz & rij_vec, xyz & fij);
    void    get_single_particle_contributions(const vector<xyz> & coords, int selected_atom, const xyz & selected_atom_coords, const xyz & boxdim, double & energy_selected, xyz & stress_selected);
    
    LJ_model(double in_sigma, double in_epsilon, double in_rcut);
    ~LJ_model();
};

LJ_model::LJ_model(double in_sigma, double in_epsilon, double in_rcut)
{
    sigma   = in_sigma;
    epsilon = in_epsilon;
    rcut    = in_rcut;
}
LJ_model::~LJ_model(){}

double LJ_model::get_eij(double rij)
{
    /* Write code to compute the Lennard Jones energy between the two atoms.
    It should account for the user-specified cutoff distance.
    The function should return the computed energy.
    */
    if(rij >= rcut){
        return 0;
    }

    double sigmaOverR = sigma/rij;

    return 4*epsilon*(pow(sigmaOverR, 12)-pow(sigmaOverR, 6));
}

void LJ_model::get_fij(double rij, const xyz & rij_vec, xyz & fij)
{
    /* Write code to compute the Lennard Jones force between the two atoms.
    It should account for the user-specified cutoff distance.
    The function should update the the force and distance vectors directly
    The function should not return anything.
    */

    if(rij >= rcut){
        fij.x = 0;
        fij.y = 0;
        fij.z = 0;
        return;
    }

    //Force is the derivative of energy with respect to R!
    //units below in J/A!!!
    fij.x = -4*epsilon*((-12*pow(sigma, 12)*pow(rij_vec.x, -13))+(6*pow(sigma, 6)*pow(rij_vec.x, -7)));
    fij.y = -4*epsilon*((-12*pow(sigma, 12)*pow(rij_vec.y, -13))+(6*pow(sigma, 6)*pow(rij_vec.y, -7)));
    fij.z = -4*epsilon*((-12*pow(sigma, 12)*pow(rij_vec.z, -13))+(6*pow(sigma, 6)*pow(rij_vec.z, -7)));

    //convert to J/m = N:
    fij.x = fij.x/A2m;
    fij.y = fij.y/A2m;
    fij.z = fij.z/A2m;

    return;
}

void LJ_model::get_single_particle_contributions(const vector<xyz> & coords, int selected_atom, const xyz & selected_atom_coords, const xyz & boxdim, double & energy_selected, xyz & stress_selected)
{
    /* Write code to determine the contributions to the total system energy, forces, and stresses due to the selected atom.
    Self interactions should not be included.*/

    energy_selected   = 0;
    stress_selected.x = 0;
    stress_selected.y = 0;
    stress_selected.z = 0;
    xyz force;
    force.x = 0;
    force.y = 0;
    force.z = 0;
    
    static double rij;
    static xyz    rij_vec;

    // Loop over all atoms. Within the loop:
    for (int i = 0; i < 0; i++)
    {
        if(i == selected_atom){
            continue;
        }
        
        // Get the scalar distance and distance vector between atoms, using MIC
        rij = get_dist(selected_atom_coords, coords[i], boxdim, rij_vec);

        // Determine pair energy, but only if interaction is within cutoff idstance
        energy_selected += get_eij(rij);

        // Determine the atom pair's contribution to the total system pressure - again, only perform 
        // if within the model cutoff - we'll use this if the move is accepted
        get_fij(rij, rij_vec, force);

        stress_selected.x += force.x * rij_vec.x * A2m;
        stress_selected.y += force.y * rij_vec.y * A2m;
        stress_selected.z += force.z * rij_vec.z * A2m;
    }
}


int main(int argc, char* argv[])
{
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Set user-defined variables (Read in from input file at some point)
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    int     seed    = stoi(argv[1]);// Seed for random number generator - read from commandline 
    double  redden  = stod(argv[2]);// Reduced density - read from commandline
    
    //////// Static variables 
    
    MTRand  mtrand(seed); // Seed the (psuedo) random number generator
    
    double  sigma   = 3.4;          // Units: Angstrom
    double  epsilon = 120;          // Units: K*k_B
    double  rcut    = 4*sigma;      // Model outer cutoff
    
    int     natoms  = 500;          // Number of atoms in the system
    string  atmtyp  = "Ar";         // Atom type (atomic symbol)
    double  molmass = 39.948;       // Mass of an particle of atmtyp (Atomic mass in the case of an argon atom)
    double  density = redden / pow(A2m,3.0) * pow(cm2m,3.0) / nav * molmass / pow(sigma,3.0);     // System density; Units: g/cm^3
    double  numden;                 // System number density; Units: atoms/Ang^3 - will be calculated later on

    double  temp    = 1.2*epsilon;  // Temperature, in K
    double  nsteps  = 5e6;          // Number of MC steps
    int     iofrq   = 2e4;          // Frequency to output statistics and trajectory
    int     nequil  = 2e6;          // Equilibration period (chemical potential and heat capacity only collected after this many steps)
    
    //////// Print information for user
    
    cout << "# Number of atoms:       " << natoms     << endl;
    cout << "# Atom type:             " << atmtyp     << endl;
    cout << "# Molar Mass (g/mol):    " << molmass    << endl;
    cout << "# Density (g/cm^3):      " << density    << endl;
    cout << "# LJ sigma (Angstrom):   " << sigma      << endl;
    cout << "# LJ epsilon/kB (K):     " << epsilon    << endl;
    cout << "# LJ cutoff (Angstrom):  " << rcut       << endl;
    cout << "# LJ cutoff (sigma):     " << rcut/sigma << endl;
    cout << "# Temperature (K):       " << temp       << endl;
    cout << "# Number MC steps:       " << nsteps     << endl;
    cout << "# Output frequency:      " << iofrq      << endl;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Initialize the system - arguments are natoms, atom type, density (g/cc), molar mass, 
    // and whether to print a simulation trajectory for object
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    system_coordinates system(natoms, atmtyp, density, molmass);
    system.generate_coords();
    
    //double check this!
    cout << "# reduc. density (rho*): " << redden << endl;
    cout << "# reduc. temp (T*):      " << (temp*kB)/epsilon << endl;

    // Intialize the model - arguments are sigma (Angstroms), epsilon (K), and outer cutoff (Angstroms)
    
    LJ_model LJ(sigma, epsilon, rcut);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Determine initial system energy and pressure
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    double  rij;                // Scalar distance between atoms
    xyz     rij_vec;            // Distance vector between atoms
    
    double  energy = 0;         // Energy for the configuration
   
    double widom_factor = 0;    // The exponential term in the Widom chemical potential expression
    double widom_trials = 0;    // Number of Widom insertion attempts
    double widom_avg    = 0;    // Average value of the wWidom factor
        
    double Eavg = 0;            // Average energy
    double Esqavg = 0;          // Square of average energy

   
    xyz     force;              // Force on a given atom
    xyz     stensor;            // System stress tensor (diagonal components only, i.e., xx, yy, zz)
    stensor.x = 0;
    stensor.y = 0;
    stensor.z = 0;
    
    for (int i=0; i<system.natoms; i++)                                                                         
    {
        for (int j=i+1; j<system.natoms; j++)
        {
            // Get the scalar distance and distance vector between atoms, using MIC
            rij = get_dist(system.coords[i], system.coords[j], system.boxdim, rij_vec);

            // Determine atom pair's contirbution to total system energy - remember to only perform the 
            // calculation if the pair distance is within the model cutoff
            energy_ij = get_eij(rij);
            energy += energy_ij;

            // Determine the atom pair's contribution to the total system pressure - again, only perform 
            // if within the model cutoff
            get_fij(dist, rij_vec, force);

            stensor.x += force.x * rij_vec.x;
            stensor.y += force.y * rij_vec.y;
            stensor.z += force.z * rij_vec.z;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Begin the simulation
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    double max_displacement = 0.5*system.boxdim.x;  // Start trial displacements at one angstrom

    int     selected_atom;
    
    double  eold_selected;  // Pre-trial-move energy of the selected atom
    double  enew_selected;  // Post-trial-move energy of the selected atom
    double  delta_energy;   // Difference between old and new (trial) energy
    xyz     sold_selected;  // Pre-trial-move stress tensor diagonal of the selected atom
    xyz     snew_selected;  // Post-trial-move stress tensor diagonal of the selected atom

    xyz     trial_displacement;    
    xyz     trial_position;
    
    int     naccepted_moves = 0;    // Number of accepted moves
    double  fraction_accepted;      // Fraction of attempted moves that have been accepted
    int     nrunningav_moves = 0;   // Move index for running averages (doesn't start until equilibration period has ended)

    double pressure;
    double Cv;

    double stat_avgE   = 0;
    double stat_avgEsq = 0;
    double stat_avgP   = 0;
    
    for (int i=0; i<nsteps; i++)
    {
        // Select a random particle. The syntax below shows how to use the random number generator. This generate a random integer between 0 and natoms-1
        
        selected_atom = int(mtrand()*system.natoms);
        
        // Determine contributions to the system's energy, force, and stress due to the selected atom 
        
        LJ.get_single_particle_contributions(system.coords, selected_atom, system.coords[selected_atom], system.boxdim, eold_selected, sold_selected);
        
        
        // Attempt to randomly move a particle - this is a multistep process
        
        // 1. Generate the trial **displacement** in x, y, and z - the particle should be able to move in positive
        // and negative directions, i.e., +/- the maximum displacement

        trial_displacement.x = (1 - (2*mtrand()))*max_displacement;
        trial_displacement.y = (1 - (2*mtrand()))*max_displacement;
        trial_displacement.z = (1 - (2*mtrand()))*max_displacement;
        
        // 2. Generate the trial **position** in x, y, and z based on the displacement
        
        trial_position.x = coords[selected_atom].x + trial_displacement.x;
        trial_position.y = coords[selected_atom].y + trial_displacement.y;
        trial_position.z = coords[selected_atom].z + trial_displacement.z;
        
        // 3. Apply PBC if the particle has moved outside the box
        
        trial_position.x -= floor(trial_position.x/system.boxdim.x)*boxdim.x;
        trial_position.y -= floor(trial_position.y/system.boxdim.y)*boxdim.y;
        trial_position.z -= floor(trial_position.z/system.boxdim.z)*boxdim.z;      
        
        // 4. Determine the energy contribution of that particle with the system **in it's trial position**
    
        LJ.get_single_particle_contributions(system.coords, selected_atom, trial_position, system.boxdim, enew_selected, snew_selected);

        if (i >= nequil) // Only do Widom tests for the equilibrated portion of the simulation
        {        
            // 5. Do a widom insertion test
            
            double ewidom;
            xyz    swidom;
            xyz    widom_position;
            
            // 5.a Generate another position for a new ghost atom
            
            widom_position.x = mtrand()*boxdim.x
            widom_position.y = mtrand()*boxdim.y
            widom_position.z = myrand()*boxdim.z
            
            
            // 5.b Calculate change in system energy due to insertion of new ghost atom 
            //is there anything we use swidom for?
            LJ.get_single_particle_contributions(system.coords, -1, widom_position, system.boxdim, ewidom, swidom);
            
            // 5.c Update the Widom factor
            
            widom_factor = exp((1/kB*temp)*(ewidom - energys))
            widom_avg   += widom_factor;
            widom_trials++;
        }
        else
        {
            // Needed to avoid a divide-by-zero during equilibration phase when these values aren't collected
            widom_avg    = 1;
            widom_trials = 1;
        }
        
        // 6. Accept or reject the move
        // If E_old is the energy of the original system and E_new is the system energy when the 
        // particle is at it's trial position, E_old - eold_selected + enew_selected = E_new
        // Therefore delta E, which = E_new - E_old is just enew_selected - eold_selected
        
        delta_energy = enew_selected - eold_selected;

        else if ( mtrand() < exp(-delta_energy/(kB*temp)) ) // Then accept
        {
            // Then the system energy has decreased **or** our random number is less than our probability to accept
            // Update the system position, energy, and stress tensor, and number of accepted moves
            system.coords[selected_atom].x = trial_position.x;
            system.coords[selected_atom].y = trial_position.y;
            system.coords[selected_atom].z = trial_position.z;

            energy += delta_energy;

            stensor.x = stensor.x - sold_selected.x + snew_selected.x;
            stensor.y = stensor.y - sold_selected.y + snew_selected.y;
            stensor.z = stensor.z - sold_selected.z + snew_selected.z;

            ++naccepted_moves;

            /*write this*/
        }

        // Update maximum diplacement to target a 50% acceptance rate
        
        fraction_accepted = float(naccepted_moves)/float(i+1);
        
        max_displacement = update_max_displacement(fraction_accepted, system.boxdim.x, max_displacement);

        // print statistics if needed - don't forget to convert the stress tensor to pressure 
        // Compute instantaneous properties
        
        //BUG - missing A2m - UNITS???
        //first term was originally numden*temp
        double pressure_kinetic = (numden/A2m)*temp*kB/natoms;
        double pressure_virial = 1.0/3.0/pow(system.boxdim.x*A2m,3.0)*(stensor.x + stensor.y + stensor.z)
        pressure = pressure_kinetic + pressure_virial;
        Cv       = 0; //calculate heat capacity
        double Cv_xs = 0;
        double Cv_ig = 0;

        if (i >= nequil) // Compute values for running averages, only using the equilibrated portion of the trajectory
        {
            stat_avgE   += energy;
            stat_avgEsq += energy*energy;        
            stat_avgP   += pressure/LJ.epsilon*pow(LJ.sigma,3.0); // Convert to reduced units! ********
            nrunningav_moves++;
            
            double avgE   = stat_avgE /float(nrunningav_moves);
            double avgEsq = stat_avgEsq / float(nrunningav_moves);
            
            double varE = avgE - avgEsq;
            Cv_xs = pow(varE, 2)/(kB*pow(temp, 2));
            Cv_ig = (3/2)*kB;
            Cv = Cv_xs + Cv_ig;
        }    
           
        if ( (i+1) % iofrq == 0)
        {
            system.write_frame(i);

            cout << "Step:  " << setw(10) << left << i;
            cout << " NAcc:  " << setw(10) << left << setprecision(3) <<  naccepted_moves;
            cout << " fAcc:  " << setw(10) << left << fixed << setprecision(3) << fraction_accepted;
            cout << " Maxd:  " << setw(10) << left << fixed << setprecision(5) << max_displacement;
            cout << " E*/N:  " << setw(10) << left << fixed << setprecision(5) << energy/natoms/LJ.epsilon;
            cout << " P*:     " << setw(10) << left << fixed << setprecision(5) << stat_avgP/float(nrunningav_moves);
            cout << " P*cold: " << setw(10) << left << fixed << setprecision(5) <<  (stat_avgP - (nrunningav_moves*pressure_kinetic))/float(nrunningav_moves);
            cout << " Mu*_xs: " << setw(10) << left << fixed << setprecision(5) << (-log(widom_factor/widom_trials)*kB*temp); 
            cout << " Cv*/N_xs:  " << setw(15) << left << fixed << setprecision(5) << Cv_xs;
            cout << " E(kJ/mol): " << setw(10) << left << fixed << setprecision(3) << energy * 0.003814; // KJ/mol per K 
            cout << " P(bar):    " << setw(10) << left << fixed << setprecision(3) << pressure * 0.003814 * 10.0e30 * 1000/(6.02*10.0e23)*1.0e-5; // KJ/mol/A^3 to bar
            cout << endl;
        }

    }
    
    stat_avgE    /= float(nrunningav_moves);
    stat_avgEsq  /= float(nrunningav_moves);
    //TODO
    Cv            =  Cv_xs/*write this, based on average energies - this should only be the dE/dT portion*/
    stat_avgE    *= 1.0/natoms/LJ.epsilon;

    cout << "# Computed average properties: " << endl;
    cout << " # E*/N:  "      << setw(10) << left << fixed << setprecision(5) << stat_avgE << endl;
    cout << " # P*:     "     << setw(10) << left << fixed << setprecision(5) << stat_avgP / float(nrunningav_moves) << endl;
    cout << " # Cv*/N_xs:   " << setw(15) << left << fixed << setprecision(5) << C/*write this - this is excess heat capacity per atom based on average energies*/ << endl;
    cout << " # Mu_xs:  "     << setw(10) << left << fixed << setprecision(5) << /*write this - this is excess chemical potential*/ << endl;       
    
}
















