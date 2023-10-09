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

    double dx = (a2.x-a1.x) - round((a2.x-a1.x)/boxdim.x) * boxdim.x;
    double dy = (a2.y-a1.y) - round((a2.y-a1.y)/boxdim.y) * boxdim.y;
    double dz = (a2.z-a1.z) - round((a2.z-a1.z)/boxdim.z) * boxdim.z;

    rij_vec.x = dx;
    rij_vec.y = dy;
    rij_vec.z = dz;

    double distance = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));

    return distance;
}

double update_max_displacement(double fraction_accepted, double boxdim, double max_displacement)
{
    /* Write code to update the maximum displacement based on current acceptance crtieria. 
    The function should return the the new maximum displacement.
    */
    if(fraction_accepted >= 0.5){
        max_displacement =  max_displacement*1.2; //if TOO MANY moves are accepted (i.e. the displacement is too short), then we want to increase max_displacement by 20%
    }

    if(fraction_accepted < 0.5){
        max_displacement = max_displacement*0.8; //if TOO FEW moves are accepted (i.e. the displacement is too large), then we want to decrease max_displacement by 20% 
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
    double      molmass;    // Molar mass of an atom of atmtyp, presumably in g/mol
    xyz         boxdim;     // Simulation box x, y, and z lengths, UNITS??
    vector<xyz> coords;     // Coordinates for all atoms in our system, UNITS??
    
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
    
    boxdim.x = pow(natoms/numden, 1.0/3.0); //ANGSTROMS
    boxdim.y = boxdim.x; //ANGSTROMS
    boxdim.z = boxdim.x; //ANGSTROMS
    
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
                
    //NOTE: atom coords is ALSO in angstroms
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
    
    double  sigma; //in ANGSTROMS
    double  epsilon; //in JOULES (or rather, kB * Kelvin)
    double  rcut; //radius - ANGSTROMS
    
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

    double sigmaOverR = sigma/rij; //unitless

    return 4*epsilon*(pow(sigmaOverR, 12)-pow(sigmaOverR, 6)); //energy returned in JOULES
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
    //power terms all have units of A^-1
    double fij_mag = -4*epsilon*((-12*pow(sigma, 12)*pow(rij, -13))+(6*pow(sigma, 6)*pow(rij, -7)));

    //convet from J/A to J/m = N
    fij_mag = fij_mag/A2m;

    fij.x = fij_mag * (rij_vec.x/rij);
    fij.y = fij_mag * (rij_vec.y/rij);
    fij.z = fij_mag * (rij_vec.z/rij);

    //FORCES RETURNED IN J/m
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
    
    //rij in ANGSTROMS
    static double rij;
    static xyz    rij_vec;

    // Loop over all atoms. Within the loop:
    for (int i = 0; i < coords.size(); i++)
    {
        if(i == selected_atom){
            continue;
        }
        
        // Get the scalar distance and distance vector between atoms, using MIC
        //rij IN ANGSTROMS
        rij = get_dist(selected_atom_coords, coords[i], boxdim, rij_vec);

        // Determine pair energy, but only if interaction is within cutoff idstance
        //eij in JOULES
        energy_selected += get_eij(rij);

        // Determine the atom pair's contribution to the total system pressure - again, only perform 
        // if within the model cutoff - we'll use this if the move is accepted
        //force in JOULES PER METER 
        get_fij(rij, rij_vec, force);

        //stress in J/m (Newtons), multiply by rij_vec (ANGSTROMS) - need to convert rij to ANGSTROMS to get JOULES
        //thus, sress is in JOULES
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
    double  epsilon = 120;          // Units: K*k_B (JOULES)
    double  rcut    = 4*sigma;      // Model outer cutoff
    
    int     natoms  = 500;          // Number of atoms in the system
    string  atmtyp  = "Ar";         // Atom type (atomic symbol)
    double  molmass = 39.948;       // Mass of an particle of atmtyp (Atomic mass in the case of an argon atom)
    double  density = redden / pow(A2m,3.0) * pow(cm2m,3.0) / nav * molmass / pow(sigma,3.0);     // System density; Units: g/cm^3
    //double  numden;                 // System number density; Units: atoms/Ang^3 - will be calculated later on

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
    
    double  rij;                // Scalar distance between atoms, in ANGSTROMS
    xyz     rij_vec;            // Distance vector between atoms, in ANGSTROMS
    
    double  energy = 0;         // Energy for the configuration, in JOULES
   
    double widom_factor = 0;    // The exponential term in the Widom chemical potential expression
    double widom_trials = 0;    // Number of Widom insertion attempts
    double widom_avg    = 0;    // Average value of the wWidom factor
        
    double Eavg = 0;            // Average energy, in JOULES
    double Esqavg = 0;          // Square of average energy, in JOULES^2

   
    xyz     force;              // Force on a given atom, given in NEWTONS (J/m)
    xyz     stensor;            // System stress tensor (diagonal components only, i.e., xx, yy, zz) IN JOULES UNTIL FINAL CONVERSION
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
            double energy_ij = LJ.get_eij(rij);
            energy += energy_ij; // energy in JOULES

            // Determine the atom pair's contribution to the total system pressure - again, only perform 
            // if within the model cutoff
            LJ.get_fij(rij, rij_vec, force); //force of initialized system in NEWTONS

            stensor.x += force.x * rij_vec.x*A2m; //stensor of initialized sytem in JOULES
            stensor.y += force.y * rij_vec.y*A2m;
            stensor.z += force.z * rij_vec.z*A2m;
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
        //displacement and position in ANGSTROMS

        trial_displacement.x = (1 - (2*mtrand()))*max_displacement;
        trial_displacement.y = (1 - (2*mtrand()))*max_displacement;
        trial_displacement.z = (1 - (2*mtrand()))*max_displacement;
        
        // 2. Generate the trial **position** in x, y, and z based on the displacement
        
        trial_position.x = system.coords[selected_atom].x + trial_displacement.x;
        trial_position.y = system.coords[selected_atom].y + trial_displacement.y;
        trial_position.z = system.coords[selected_atom].z + trial_displacement.z;
        
        // 3. Apply PBC if the particle has moved outside the box using MIC:
        
        trial_position.x -= floor(trial_position.x/system.boxdim.x)*system.boxdim.x;
        trial_position.y -= floor(trial_position.y/system.boxdim.y)*system.boxdim.y;
        trial_position.z -= floor(trial_position.z/system.boxdim.z)*system.boxdim.z;      
        
        // 4. Determine the energy contribution of that particle with the system **in it's trial position**
        //energy is in JOULES, stress is ALSO in JOULES
        LJ.get_single_particle_contributions(system.coords, selected_atom, trial_position, system.boxdim, enew_selected, snew_selected);

        if (i >= nequil) // Only do Widom tests for the equilibrated portion of the simulation
        {        
            // 5. Do a widom insertion test
            
            double ewidom;
            xyz    swidom;
            xyz    widom_position;
            
            // 5.a Generate another position for a new ghost atom
            
            widom_position.x = mtrand()*system.boxdim.x;
            widom_position.y = mtrand()*system.boxdim.y;
            widom_position.z = mtrand()*system.boxdim.z;
            
            
            // 5.b Calculate change in system energy due to insertion of new ghost atom 
            //is there anything we use swidom for?
            LJ.get_single_particle_contributions(system.coords, -1, widom_position, system.boxdim, ewidom, swidom);
            
            // 5.c Update the Widom factor
            //check this! (deltaE portion)
            
            widom_factor = exp((1/kB*temp)*(ewidom)); //unit check: ewidom is in JOULES, beta is in 1/(J/K * K), thus in 1/J, so the exp factor is unitless (yes!)
            widom_avg   += widom_factor; //widom factor, and widom average, are unitless
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

        if ( mtrand() < exp(-delta_energy/(kB*temp)) ) // Then accept
        {
            // Then the system energy has decreased **or** our random number is less than our probability to accept
            // Update the system position, energy, and stress tensor, and number of accepted moves
            //update selected atom's location
            system.coords[selected_atom].x = trial_position.x;
            system.coords[selected_atom].y = trial_position.y;
            system.coords[selected_atom].z = trial_position.z;

            //add the change in energy caused by this particle perturbation to the total energy
            energy += delta_energy;

            //update stress tensor - subtract the old stress and add the new stress!
            stensor.x = stensor.x - sold_selected.x + snew_selected.x;
            stensor.y = stensor.y - sold_selected.y + snew_selected.y;
            stensor.z = stensor.z - sold_selected.z + snew_selected.z;

            //update the number of updated moves!    
            ++naccepted_moves;

            /*write this*/
        }

        // Update maximum diplacement to target a 50% acceptance rate
        
        //calculated the number of accepted moves - number of accepted moves divided by the number of loops (i+1 because C++ is zero indexed)
        fraction_accepted = float(naccepted_moves)/float(i+1);
        
        //update the max displacement to target a 50% acceptance rate 
        max_displacement = update_max_displacement(fraction_accepted, system.boxdim.x, max_displacement);

        // print statistics if needed - don't forget to convert the stress tensor to pressure 
        // Compute instantaneous properties
        
        //BUG - missing A2m - UNITS???
        //first term was originally numden*temp
        double pressure_kinetic = (system.numden/A2m)*temp*kB/natoms; //UNIT CHECK - numdem is in atoms / A^3 - divide by A2m to get atoms/m^3, multiply by K * J/K, so multiply by ENERGY (J * atoms / m^3), then divide by the number of atoms to get J/ m^3
        double pressure_virial = 1.0/3.0/pow(system.boxdim.x*A2m,3.0)*(stensor.x + stensor.y + stensor.z); //UNIT CHEKC - stensors are all in JOULES - divide by volume of the box, which is the box.dim cubed, and is in ANGSTROMS, so convert to meters and cube - units in J/m^3
        pressure = pressure_kinetic + pressure_virial; //PRESSURE IN J/m^3
        Cv       = 0; //calculate heat capacity
        double Cv_xs = 0;
        double Cv_ig = 0;

        if (i >= nequil) // Compute values for running averages, only using the equilibrated portion of the trajectory
        {
            stat_avgE   += energy; //energy in JOULES
            stat_avgEsq += energy*energy;        
            stat_avgP   += pressure/LJ.epsilon*pow(LJ.sigma*A2m,3.0); // Convert to reduced units! ******** pressure in J/m^3, divide by joules, multiply by sigma^3** sigma is in angstroms, so must convert to meters!!
            nrunningav_moves++;
            
            double avgE   = stat_avgE /float(nrunningav_moves);
            double avgEsq = stat_avgEsq / float(nrunningav_moves);
            
            double varE = avgEsq - pow(avgE, 2); //varE is in J^2
            Cv_xs = varE/(kB*pow(temp, 2)); //J^2, divided by J/K * K^2, leaves Cv units to J/k (yes!)
            //Cv_ig = (3/2)*kB;
        }    
           
        if ( (i+1) % iofrq == 0)
        {
            system.write_frame(i);

            cout << "Step:  " << setw(10) << left << i;
            cout << " NAcc:  " << setw(10) << left << setprecision(3) <<  naccepted_moves;
            cout << " fAcc:  " << setw(10) << left << fixed << setprecision(3) << fraction_accepted;
            cout << " Maxd:  " << setw(10) << left << fixed << setprecision(5) << max_displacement; 
            cout << " E*/N:  " << setw(10) << left << fixed << setprecision(5) << energy/natoms/LJ.epsilon; //energy per atom in reduced units - epsilon in JOULES
            cout << " P*:     " << setw(10) << left << fixed << setprecision(5) << stat_avgP/float(nrunningav_moves); //total average pressure, in reduced units 
            cout << " P*cold: " << setw(10) << left << fixed << setprecision(5) <<  (stat_avgP/float(nrunningav_moves) - pressure_kinetic/LJ.epsilon*pow(LJ.sigma*A2m, 3.0));
            cout << " Mu*_xs: " << setw(10) << left << fixed << setprecision(5) << (-log(widom_avg/widom_trials)*kB*temp); //excess chemical potential
            cout << " Cv*/N_xs:  " << setw(15) << left << fixed << setprecision(5) << Cv_xs/natoms; //heat capacity per atom, in J/(K*atom)
            cout << " E(kJ/mol): " << setw(10) << left << fixed << setprecision(3) << energy/natoms*nav/1000; // KJ/mol per K 
            cout << " P(bar):    " << setw(10) << left << fixed << setprecision(3) << pressure / natoms / pow(10, 5); // KJ/mol/A^3 to bar
            cout << endl;
        }

    }
    
    stat_avgE    /= float(nrunningav_moves);
    stat_avgEsq  /= float(nrunningav_moves);
    //TODO - calculate the average heat capacity based on average energy???
    Cv            =  (stat_avgEsq - pow(stat_avgE, 2))/(kB*pow(temp, 2));/*write this, based on average energies - this should only be the dE/dT portion*/
    stat_avgE    *= 1.0/natoms/LJ.epsilon;

    cout << "# Computed average properties: " << endl;
    cout << " # E*/N:  "      << setw(10) << left << fixed << setprecision(5) << stat_avgE << endl;
    cout << " # P*:     "     << setw(10) << left << fixed << setprecision(5) << stat_avgP / float(nrunningav_moves) << endl;
    cout << " # Cv*/N_xs:   " << setw(15) << left << fixed << setprecision(5) << Cv/natoms << endl; //excess heat capacity per atom
    cout << " # Mu_xs:  "     << setw(10) << left << fixed << setprecision(5) << -log(widom_avg/widom_trials)*kB*temp << endl; //excess chemical potential    
    
}
















