// perp-parallel-joel.cpp
//
// Simulates Hbar trajectories in the ALPHA trap
// -Record specific trajectories (x(t), v(t) for all timesteps)
// -Record E_{\parallel}(t) whenever z(t) crosses 0
// -Compute largest Lyapunov exponents (and save them whenever necessary)
// -Records threshold times
//
// Mike Zhong
// 7 June, 2017



//standard header files
#include <string.h> //string manipulation
#include <stdio.h>   // file I/O
//#include <math.h>    // alternate math functions
#include "/usr/include/math.h" // must use this for running on school's cluster
#include <time.h>    // time functions
#include <stdlib.h>  // added to handle calls to getenv()
#include <iostream> // FOR PRINTING OUTPUT FILES
#include <fstream> // to write output files
#include <vector> // arrays of constant
#include <ctime> // for getting current time, to write
#include <sstream> // for string manipulations


//---------------------------------- MPI library/instructions for parallelization
//standard parallelization header files
// MPI UNCOMMENT
// #include <mpi.h>
//----------------------------------


// MAIN PARAMETERS: EDIT HERE
int numtraj = 100; //number of trajectories to simulate
double tfin = 100.0; // final time of each trajectory
double dtim = 3.5e-5; //dtim is the step size of the symplectic stepper.
                      //3.5e-7 is the usual value for the simulations.

char* potential_option; // "oct", "quad", or "ho"
char* test_type; // "lle", "thresh", or "crossings"

// parameters for Largest Lyapunov Exponent
double dx0 = 1.0e-8; // size of perturbation for calculating lyapunov exponent
// we want v_scale to be s.t. at velocity s.t. kinetic energy is 0.54 Kelvins, the conversion would be the length of trap
// 0.022275m = 22.275 cm
// (m v_{max} ^2)/2 = 0.54 K * k_B = 7.45550201 × 10^-24 joules
// m_{prot} = 1.672621777E-27 kg; v_{max} = sqrt( 2 * 7.45550201 × 10^-24 joules / 1.672621777E-27 kg) = 94.4179507 m / s
// conversion to v_{max}: v_scale * 94.4179507 m / s = 0.022275m
// v_scale = 0.000235919122
double v_scale = 0.000235919122; // scaling for velocity in metric len((x,v)) = sqrt(x^2 + (v_scale*v)^2)
double init_pert[6] = {dx0, 0, 0, 0, 0, 0}; // initial perturbation


// parameter for threshold
double thresh_frac = 0.2;
double measure_time = 0.0; // start measurement at this time

// parameters used in main and in subroutines
// [SI units + Kelvin wherever specified]
const double pi = double(2)*asin(double(1)),
kboltz = 1.3806488E-23,
mprot = 1.672621777E-27,
sqrt12th = double(1)/sqrt(double(12)),
magmom0 = 9.27400968E-24;

// Other constants related to magnetic field calculations.
// treat the mirror coils as 2 solenoids that start/end at zsol_0/zsol_f
// compared to the center position
const double zsol10 = -14.292e-3, zsol1f = 14.292e-3, lam1 = 0.90190,
             zsol20 = -14.859e-3, zsol2f = 14.859e-3, lam2 = 0.90266,
             radis = 4.5238e-2;
const double z1 = 0.12946, imu01 = -2.1632e4/1100, dz1 = 16.449e-3;

// magnetic field things
double radtrp, radtrp2; // radius of trap
double sm_rad, sm_rad2; // radius of smaller neck (near ends)
double md_rad, md_rad2; // radius medium part of trap
double end1, end2, end3, end4, end5; // variables to do with annihilation with walls

// random number generator parameters
const long int ntab = 16;
const long int ndiv2 = 1+(2147483563-1)/ntab; // 134217723.625; 2147483563 ~ 2^31; 2147483648 = 2^31
long int iv1[ntab];

double zmid; // middle of trap

// Global variables related to magnetic field calculations
double bfmin;          // minimum B-field value, to compute energy
double bfedge;         // B-field at edge of trap
double bunif;          // uniform B-field
double bsol1, bsol2;   // B-field from mirror coils on ends
double co1, co2, co3, co4, co5, co6; // octopole coefficients
double imu1, ioct;     // more octopole values
double imu2;           // quadropole values

// coefficients for the mirror field calculation
double zmirr1, zmirr2; //zmirr1 and zmirr2 are the z-positions of the mirror coils
double bcentsol1, bcentsol2, imu0a1, imu0a2, radis2, zpos11, zpos12, zpos21, zpos22;

// ---------------------- FUNC DECLARATIONS! ----------------------
// initializes the random number generator
void InitRan(long int& i1, long int& i2);

// computes random #
double SimRan1(long int& i1, long int& i2);

// computes the complementary error function using a fit good to ~1E-7
//  from Numerical Recipes
double Erfcc(double x);

// computes the B-field from a solenoid
void Bsole(double radi, double zpos0, double zposf, double bsole, double r[], double bf[]);

// computes the B-field from a solenoid uses 2 simplifications for B_sole
//   that vastly speeds up the calculation but can be problematic near
//   the coils
void Bsolesimp(double r[], double bf[]);

// computes the B-field from a circle loop
void Bcirc(double radi, double zpos, double bcent, double r[], double bf[]);

// computes the actual B_oct field using a fit to the ALPHA fields inside
//   the trap regions; it gives accuracy of 1% for z>0 but I did not
//   know how to do the leads
void Boct(double r[], double bf[]);

// computes the actual B_quad field
void Bquad(double r[], double bf[]);

// computes the B_oct field without edges
void Boct_edgeless(double r[], double bf[]);

// computes net B-field
// includes the 2 mirror coils, the uniform field, and the octupole field
void Bfield(double r[], double bf[]); // changed name from Bfield_mod

// computes MAGNETIC force on Hbar. It uses the currents calculated in UpdateCurrents_mod()
void ForceB(double r[], double frc[]);

// computes parameters needed to compute the potential and force
void InitialSetup();

// generates initial conditions for Hbar
void GenInitCond(long int& i1, long int& i2, double xv[]);

// symplectic stepper
// takes xv[t] -> xv[t + dtim]
void Symplec(double t, double dtim, double xv[]);

// compute the time-derivatives of the variables
// used in initialization, by pushing x(t) back by half a time-step
void Dydt(double t, double xv[], double dxvdt[]);

// calculates energy of Hbar given its position and velocity
double Energy(double xv[]);

// computes PE of Pbar in external fields
double PotEnr(double x, double y, double z);

// computes DOT PRODUCT
double Dot(double x[], double y[]);

// computes CROSS PRODUCT
double * Cross(double x[], double y[]);

// decides whether to continue iteration loop
// depends on whether particle is still within trap, and whether t < tfin
bool ContinueLoop(double tim, double tfin, double xv[3]);

// returns current date and time as string
std::string DateString();

// writes to individual files each trajectory
void WriteMetaParams(std::string name, int numtraj, double tfin, double dtim, int seed);

// writes to output file the crossings information
void WriteCrossings(std::string name, std::vector<double> cross_t, std::vector<double> cross_Tz);

// computes the average of a file -- returns the average T_z value
long double average(std::vector<double> cross_Tz);

// computes square average, for variance calculation
long double variance(std::vector<double> cross_Tz, long double average);



//  ---------------------------------- MAIN  ----------------------------------
int main(int argc, char *argv[]) {

  int inode, nproc;   // MPI variables
  inode = 1; // overwritten by MPI
  int seed; // random algorithm coefficients!!
  long int i1, i2;
  std::string init_name;

  int itraj; // iterating over number of trajectories to be taken

  double tim; // time of simulation
  double timetaken; // computation time to do simulation
  double xv[6], xv0[6], dxvdt[6], energy0; // phase-space variables, and init energy

  // perturbation variables for LLE calculation
  double xv_pert[6];
  double pert[6];
  double dx_;

  time_t todstr, todfin, todstr2;
	time(&todstr);

  /*
  //----------------------------------
  // MPI UNCOMMENT

    //initialize mpi variables
    MPI_Init(&argc, &argv);
    //find out how many processors
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    //find out which node
    MPI_Comm_rank(MPI_COMM_WORLD, &inode);
    printf("Num. processors: %d, Num. node: %d\n", nproc, inode);

  //----------------------------------
  */

  //seed = atoi(getenv("SEED"));             // Get "seed" from the environment.

  if (argc != 4) {       // Data file is passed from the command line.
    printf("Wrong number of options!\nUsage: ./perp-parallel-joel potential_option test_type seed_number\nTry harder next time!\n");
    return 1;
  } else {
    potential_option = argv[1];
    test_type = argv[2];
    seed = atoi(argv[3]);

    if (strcmp(potential_option, "oct") == 0) {
      printf("Normal Octupole\n");
    } else if (strcmp(potential_option, "quad") == 0) {
      printf("Quadrupole\n");
    } else if (strcmp(potential_option, "ho") == 0) {
      printf("Harmonic Oscillator!\n");
    } else {
      printf("That is not a good potential_option.\n");
      printf("Options are: oct, quad, or ho.\n");
      printf("Better luck next time!\n");
      return 2;
    }

    if (strcmp(test_type, "all") == 0) {
      printf("To be printed out: all information (including LLE)\n");
      printf("t,x,y,z,perp_x,perp_y,perp_z,perp_vx,perp_vy,perp_vz,dlog|perp|/dt");
    } else if (strcmp(test_type, "lle") == 0) {
      printf("Largest Lyapunov Exponent\n");
    } else if (strcmp(test_type, "cross") == 0) {
      printf("z = 0 Crossings\n");
    } else if (strcmp(test_type, "thresh") == 0) {
      printf("Threshold Measurements!\n");
    } else {
      printf("That is not a good test type.\n Options are: lle, cross, or thresh");
      return 3;
    }

    if (seed < 1) {
      printf("Seed must be positive\nTry harder next time!\n");
      return 4;
    }
  }

  init_name = DateString() + "_" + potential_option + "_" + test_type;

  // MPI UNCOMMENT
  // seed += inode;         // Generate different seeds for the parallel processes.
  // sleep(inode*10);       // Delay to test that the MPI execution is not messed up.

  if(seed < 0) seed = (-1)*seed; seed ++;    // Make sure that seed is positive.

  // initialize the random number generator
  i1 = 5+2836*(seed+115); i2 = 123+12211*(seed+115);
  InitRan(i1,i2);

  // Setup the information for the electrodes and the magnetic fields.
  InitialSetup();


  // ------------------- START THE SIMULATION!! ------------------------

  if (inode == 1) {
    WriteMetaParams(init_name, numtraj, tfin, dtim, seed);
  }

  // initializes params file, in which each trajectories' parameters are written
  std::ofstream params;
  params.open(("params_" + init_name + ".csv").c_str(), std::ios_base::app);
  params.close();


  // START TRAJECTORY
  time(&todstr2);
  itraj = 0;

  while(itraj < numtraj) {  //printf("\nn = %d ", itraj);
    int cnt = 0;
    int lyapunov_n = 0;

    // initialize the time, starting energy and potentials at the electrodes.
    tim = -1.0 + 0.1*(2.0*SimRan1(i1,i2) - 1.0);    // start time between -1.1 and -0.9 seconds
    GenInitCond(i1, i2, xv);		// Generate random initial conditions for this iteration of the loop.

    // save the initial conditions for output of trajectories that survive.
    for (int j = 0; j < 6; j++) {
      xv0[j] = xv[j];
    }
    energy0 = Energy(xv);
    // printf("\nE0=%15.8E ", energy0 / kboltz);

    // shift the position backward in time by 1/2 of time step for Leapfrog integration

    Dydt(tim,xv,dxvdt);

    // taylor expansion, x(t + dt) ~= x(t) + dt*x'(t) + 0.5*(dt)^2*x''(t)
    double dtim_init = -0.5*dtim;
    for (int j = 0; j < 3; j++) {
      xv[j] += dxvdt[j]*dtim_init + 0.5*(dtim_init*dtim_init)*dxvdt[j+3];
    }
    // now, x is at tim - dtim/2, v is at tim

    // initialize Lyapunov exponent algorithm
    for (int j = 0; j < 6; j++) {
      xv_pert[j] = xv[j] + init_pert[j];
      pert[j] = init_pert[j];
    }

    if (strcmp(test_type, "all") == 0) {
      printf("\n%15.8E,", tim);
      for (int j = 0; j < 3; j++) {
        printf("%15.8E,", xv[j]); // only need to print x, y, z; vx, vy, vz can be determined from data
      }
      for (int j = 0; j < 6; j++) {
        printf("%15.8E,", pert[j]); // only need to print x, y, z; vx, vy, vz can be determined from data
      }
      printf("0,"); // start with 0 in d(log pert/dt))
    }

    // initialize z = 0 crossings
    double xv_prev[6], xv_prev2[6];
    double last_z, curr_z;
    double T_z; // at z ~= 0, axial potential energy should be nearly zero
    double vz_m1, vz_0, vz_p1;
    double vz_cross;
    double s;
    double xv_cross[6];
    int cross_n;

    std::vector<double> cross_t;
    std::vector<double> cross_Tz;
    std::vector<double> cross_x;
    std::vector<double> cross_y;
    std::vector<double> cross_vx;
    std::vector<double> cross_vy;


    double lyapunov_sum = 0;


    // threshold parameters
    double Tz_init;
    double t_init;
    bool started_thresh = false;
    bool finished_thresh = false;

    cross_n = 0;
    while ((ContinueLoop(tim, tfin, xv)) && (!(finished_thresh))) { // finished_thresh only set in threshold test
      // should step x(t - dt/2) -> x(t + dt/2); v(t) -> v(t+dt);
      Symplec(tim, dtim, xv); // time-evolves xv by one step

      //increment time and total number of time steps
      tim += dtim;
      cnt++;

      if (strcmp(test_type, "all") == 0) {
        printf("\n%15.8E,", tim);
        for (int j = 0; j < 3; j++) {
          printf("%15.8E,", xv[j]); // only need to print x, y, z; vx, vy, vz can be determined from data
        }
      }
          // printf("E0=%15.8E ", energy0 / kboltz);


      if ((strcmp(test_type, "lle") == 0) || (strcmp(test_type, "all") == 0)) {// Lyapunov Algorithm
        Symplec(tim, dtim, xv_pert); // time-evolves xv_pert by one step

        // computing perturbation size, dx_
        for (int j = 0; j < 6; j++) {
          pert[j] = xv_pert[j] - xv[j];
          if (strcmp(test_type, "all") == 0)
            printf("%15.8E,", pert[j]); // print out the new perturbation
        }
        dx_ = 0;
        for (int j = 0; j < 3; j++) {
          dx_ += pert[j]*pert[j];
        }
        for (int j = 3; j < 6; j++) {
          dx_ += (v_scale*pert[j])*(v_scale*pert[j]);
        }
        dx_ = sqrt(dx_);

        if (strcmp(test_type, "all") == 0)
          printf("%15.8E,", log(dx_ /dx0) / dtim); // print out the log perturbation growth rate

        // adds to lyapunov_sum the log perturbation growth
        lyapunov_sum += log(dx_ /dx0);

        // rescales perturbation to be magnitude dx0
        for (int j = 0; j < 6; j++) {
          xv_pert[j] = xv[j] + pert[j]*(dx0 / dx_);
        }
        /*
        if (cnt % 100000 == 0) {
          printf("(%15.8E, ", tim);
          printf("%15.8E), ", lyapunov_sum / (cnt * dtim));
        }*/
        /*
        printf("dx=%15.8E\n", (dx_ / dx0));
        printf("lyapunov_sum = %15.8E\n", lyapunov_sum);
        */

      } else if ((strcmp(test_type, "cross") == 0) || (strcmp(test_type, "thresh") == 0)) {
        curr_z = xv[2];
        last_z = xv_prev[2];

        if (((last_z < 0) && (curr_z >= 0)) || ((last_z >= 0) && (curr_z < 0))) {
          /* s is linear extr. between (t - dt/2) and (t + dt/2) for z(t + s*dt) = 0
            s (between -1/2 and 1/2) that satisfies:
            0 = (1/2 - s)*(last_z) + (s + 1/2)*curr_z = 0.5(last_z + curr_z) - s(last_z - curr_z)
            => s = 0.5*(curr_z + last_z) / (last_z - curr_z)
          */
          s = 0.5*(curr_z + last_z) / (last_z - curr_z);

          //printf("s=%15.8E, \n", s);

          //printf("z_cross=%15.8E, \n", (0.5 - s) * xv_prev[2] + (s + 0.5)*xv[2]);

          // xv contains x(t + dt/2), v(t + dt)
          // xv_prev contains x(t - dt/2), v(t)
          // xv_prev contains x(t - 3*dt/2), v(t - dt)

          // we care only about the vz's:
          // vz[t - dt] = vz[-1]; vz[t] = vz[0]; vz[t + dt] = vz[1]
          // quadratic interpolation
          vz_m1 = xv_prev2[5];
          vz_0 = xv_prev[5];
          vz_p1 = xv[5];

          // vz(s) = vz_m1*(s)(s-1)/(-1)(-1-1) + vz_0*(s+1)(s-1)/(0+1)(0-1) + vz_p1*(s+1)(s)/(1+1)(1)
          vz_cross = vz_m1*s*(s-1)/2 + vz_0*(s-1)*(s+1)/(-1) + vz_p1*s*(s+1)/2;

          for (int i = 0; i < 3; i++) {
            xv_cross[i] = (0.5 - s)*xv_prev[i] + (s + 0.5)*xv[i];
          }
          for (int i = 3; i < 6; i++) {
            xv_cross[i] = xv_prev2[i]*s*(s-1)/2 + xv_prev[i]*(s-1)*(s+1)/(-1) +
                          xv[i]*s*(s+1)/2;
          }

          T_z = 0.5*mprot*(vz_cross*vz_cross);

          if (strcmp(test_type, "cross") == 0) { // only for cross test
            cross_t.push_back(tim + dtim*s); // more accurate ...
            cross_Tz.push_back(T_z / kboltz);
            cross_x.push_back(xv_cross[0]);
            cross_y.push_back(xv_cross[1]);
            cross_vx.push_back(xv_cross[3]);
            cross_vy.push_back(xv_cross[4]);
          }

          else if (strcmp(test_type, "thresh") == 0) { // only for threshold test_type
            if (!(started_thresh) && ((tim + s*dtim) > measure_time)) { // take initial measurement
              started_thresh = true;
              t_init = tim + s*dtim;
              Tz_init = T_z;
              bool started_thresh = false;
              bool finished_thresh = false;
            }
            else if (started_thresh) {
              if ((T_z - Tz_init > energy0*thresh_frac) || (T_z - Tz_init < -energy0*thresh_frac)) {
                finished_thresh = true; // finished sequence
                printf("\n(t_init,Tz_init,t_fin,Tz_fin)");
                printf("\n(%8.8E,", t_init);
                printf("%8.8E,", Tz_init);
                printf("%8.8E,", (tim + s*dtim));
                printf("%8.8E)", T_z);
              }
            }
          }

          // increases number of crosses
          cross_n++;
        }

        // updating xprev
        for (int i = 0; i < 6; i++) {
          xv_prev2[i] = xv_prev[i];
          xv_prev[i] = xv[i];
        }
      }
    }
    // jot down parameter file
    params.open(("params_" + init_name + ".csv").c_str(), std::ios_base::app);
    params << seed << "," << itraj << "," << init_name << "," << energy0 / kboltz << ","
    << xv0[0] << "," << xv0[1] << "," << xv0[2] << "," << xv0[3] << "," << xv0[4] << "," << xv0[5] << "\n";
    params.close();

    if ((strcmp(test_type, "thresh") == 0) && (!(finished_thresh))) {
      printf("\n(t_init,Tz_init,t_fin,Tz_fin)");
      printf("\n(%8.8E,", t_init);
      printf("%8.8E,", Tz_init);
      printf("nan,nan)");
    }

    printf("\n(%15.8E,", energy0 / kboltz);
    printf("%15.8E", tim);

    if ((strcmp(test_type, "lle") == 0) || (strcmp(test_type, "all") == 0))
      printf(",%15.8E),\n", lyapunov_sum / (cnt*dtim));
    else
      printf(")");
    //printf("calculated lyapunov exponent=%15.8E", lyapunov_sum / (dtim * cnt));
    // shift position forward in time by 1/2.
    // to be honest ... final position doesn't matter ... does it?
    /*Dydt(tim,xv,dxvdt);
    double dtim_fin = +0.5*dtim;
    for (int j = 0; j < 3; j++) {
        xv[j] += (dxvdt[j]*dtim_fin) + 0.5*(dxvdt[j+3]*dtim_fin*dtim_fin);
    } */


    if ((strcmp(test_type, "cross") == 0) && (cross_n > 0)) { // as long as there is at least one crossing
      std::string name;
      std::stringstream name_sstream;
      name_sstream << "seed";
      name_sstream << seed;
      name_sstream << "_";
      name_sstream << DateString();
      name_sstream >> name;

      // reference: params << "seed,itraj,name,energy,total_time,cross_n,corr_int\n";
      //params.open(("params_" + init_name + ".csv").c_str(), std::ios_base::app);
      //params << seed << "," << itraj << "," << name << "," << energy0 / kboltz << "," << (tim  - cross_t[0]) << "," << cross_n << "," << corr_int << "," << avg << "," << var <<"\n";
      //params.close();
      WriteCrossings(name, cross_t, cross_Tz);
    }



    itraj++;
  }

  time(&todfin);
  printf("Total simtime is (s): %13.6E\n", difftime(todfin,todstr2));

  /*
  //------------------------
  // MPI UNCOMMENT
   //barrier stops the computation until all nodes reach this point
    MPI_Barrier(MPI_COMM_WORLD) ;

    //finalize the mpi variables before the end of the program
    MPI_Finalize();
  //---------------
  */
  return 1;
} // END OF MAIN()


//--------- takes the dot product
double Dot(double a[], double b[]) {
  double dot = 0.0;
  for (int j = 0; j < 3; j++) {
    dot += a[j]*b[j];
  }
  return dot;
}

//--------- takes the cross product
void Cross(double a[], double b[], double axb[]) {
  axb[0] = a[1]*b[2] - a[2]*b[1];
  axb[1] = a[2]*b[0] - a[0]*b[2];
  axb[2] = a[0]*b[1] - a[1]*b[0];
}

//--------- calculates the error function!!
double Erfcc(double x) {
	double t, z, ans;
	z = fabs(x);
	t = 1.0 / (1.0 + 0.5*z);
	ans = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		    t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		    t*(-0.82215223+t*0.17087277)))))))));
  if (x >= 0) {
    return ans;
  } else {
    return (2.0 - ans);
  }
}

//--------- function to initiate the random number generator
void InitRan(long int& i1, long int& i2) {
	int iran;
	for (iran = 0 ; iran < 152 ; iran++) {
  	i1 = 16807*(i1%127773) - 2836*(i1/127773);
    if (i1 <= 0)
      i1 += 2147483647;
  	i2 = 40014*(i2%53668) - 12211*(i2/53668);
    if (i2 <= 0)
      i2+= 2147483563 ;
	}
	for (iran = 0 ; iran < ntab ; iran++) {
  	i1 = 16807*(i1%127773) - 2836*(i1/127773);
    if (i1 <= 0)
      i1 += 2147483647;
  	i2 = 40014*(i2%53668) - 12211*(i2/53668);
    if (i2 <= 0)
      i2+= 2147483563 ;
  	iv1[iran] = i1; // RANDOMIZES (somehow) iv1, between 1 and 2147483563
	}
}

// -------- Function that initializes the parameters needed for the calculations of the magnetic and electric fields.
void InitialSetup() {
	int j, imid; //i, k;
	double zm1, zm2;
	double bf[3], rp[3];

	//GEOMETRIES of the trap
  radtrp = 0.022275;   // Trap radius
  radtrp2 = radtrp*radtrp;
  sm_rad = 0.014800;   // Trap radius (small electrode radius)
  sm_rad2 = sm_rad*sm_rad;
  md_rad = 0.016000;
  md_rad2 = md_rad*md_rad;

  zmid = 0;

  end1 = -0.183385; // end of trap
  end2 = -0.1442752; // between small trap and medium trap
  end3 = -0.1392752; // between medium trap and big radius
  end4 = 0.1392752; // between big trap and small radius
  end5 = 0.195095; // end of trap

	// --- Magnetic field calculation data, some of it was previously read from hbar_qdif_lin2010.datt
  bunif = 1.0;        // Uniform magnetic field strength (T)
	ioct = 885.7; //890.0; 	    // Octupole current (A); produces a 1.54T field at the wall and z=0. (890.0A in Francis' simulations)
	bsol1 = 626.5; //608.0;      // Mirror 1 current (A); produces a maximum on-axis field of 1T. (608.0A in Francis' simulations)
	bsol2 = 626.5; //608.0;      // Mirror 2 current (A); produces a maximum on-axis field of 1T. (608.0A in Francis' simulations)

  zm1 = -0.13701;     // Position of mirror coil 1 relative to the middle of the trap
	zm2 = 0.13701;      // Position of mirror coil 2 relative to the middle of the trap

	//coefficients for the octupole field calculation
	imu1 = imu01*ioct;
	co1 = 2/(dz1*sqrt(pi)); co2 = co1*2/dz1; co3 = co2*2/dz1;
	co4 = co3*2/dz1; co5 = co4*2/dz1; co6 = co5*2/dz1;


  //coefficients for the quadrupole field calculation
  imu2 = 1.54 / radtrp;

	//coefficients for the mirror field calculation
	zmirr1 = zmid+zm1; zmirr2 = zmid+zm2;    //compute the actual positions of the mirror coils
	bcentsol1 = (1.1987/750)*bsol1*sqrt(1.0+0.25*(zsol1f-zsol10)*(zsol1f-zsol10)/(radis*radis));
	bcentsol2 = (1.1929/750)*bsol2*sqrt(1.0+0.25*(zsol2f-zsol20)*(zsol2f-zsol20)/(radis*radis));
	imu0a1=2.0*radis*bcentsol1*0.25*radis*radis;
	imu0a2=2.0*radis*bcentsol2*0.25*radis*radis;
	radis2 = radis*radis;
	//mirror 1
	//compute the B-field at the center of 1 loop that gives B_sole at the center of the solenoid the position of mirror1, loop 1
	zpos11 = zmirr1+ (zsol10+zsol1f)/2 - sqrt12th*(zsol1f-zsol10);
	zpos12 = zmirr1+ (zsol10+zsol1f)/2 + sqrt12th*(zsol1f-zsol10);   	//the position of mirror1, loop 2

	//mirror 2
	//compute the B-field at the center of 1 loop that gives B_sole at the center of the solenoid the position of mirror2, loop 1
	zpos21 = zmirr2 + (zsol20+zsol2f)/2 - sqrt12th*(zsol2f-zsol20);
	zpos22 = zmirr2 + (zsol20+zsol2f)/2 + sqrt12th*(zsol2f-zsol20); 	//the position of mirror2, loop 2

	//compute the parameter "bfmin" that is used to calculate the energy in Energy_mod().
	rp[0] = 0.0; rp[1] = 0.0; rp[2] = zmid;
	Bfield(rp,bf);			// The B field is calculated using parameters for t=0.
	bfmin = sqrt(bf[0]*bf[0] + bf[1]*bf[1] + bf[2]*bf[2]);

	//compute the parameter "bfedge" that is used in the Random Initial Condition generator, GenInitCond().
	rp[0] = 0.0; rp[1] = radtrp; rp[2] = zmid;
	Bfield(rp,bf);			// The B field is calculated using parameters for t=0.
	bfedge = sqrt(bf[0]*bf[0] + bf[1]*bf[1] + bf[2]*bf[2]);   // <--- "bfedge" is needed in the calculations below
}


//--------- Generate initial conditions for Hbar
//          Taken from the loop over the total number of trajectories to simulate.
void GenInitCond(long int& i1, long int& i2, double xv[]) {
  double dum, dum2, theta, vexp, sigma;
	double Tcut, Tdist;

  //Francis initial position generator
	dum = 0.8e-3*pow(SimRan1(i1,i2),(1./3.)); // 0.0008 * R ^ (1/3) (R between 0 and 1)
	dum2 = 2.0*SimRan1(i1,i2) - 1.0; // Random between -1 and 1
  theta = 2.0*pi*SimRan1(i1,i2); // Random between 0 and 2pi

  // CONSTRAINED on ellipse defined by: r^2 + (z/10)^2 = 1, r:=radial, z:=axial
	xv[0] = dum*sqrt(1.0-dum2*dum2)*cos(theta);
  xv[1] = dum*sqrt(1.0-dum2*dum2)*sin(theta);
	xv[2] = zmid + 10.0*dum*dum2;

  // Calculation of velocities following Francis's recipe (sent by email):
	Tdist = 50.0;        // Temperature of the Hbars when they are made (~50 K)
	Tcut = 0.75;          // Minimum energy at which we know that Hbar cannot be confined in the trap (K)
  // fine tune Tcut?
	vexp = sqrt(2.0*kboltz*Tdist/mprot);
  sigma = sqrt(kboltz*Tdist / mprot);
	dum = Tcut+1.0 ;
	while(dum > kboltz*Tcut) // while K.E. of Hbar is greater than Tcut:
	{ // see Box-Muller Transformation
    // (http://mathworld.wolfram.com/Box-MullerTransformation.html)
		theta = 2.0*pi*SimRan1(i1,i2);

    // dum = vexp*sqrt(-log(SimRan1(i1,i2))) ;
    // this should give same result: (sigma = vexp/sqrt(2))
    dum = sigma*sqrt(-2.0*log(SimRan1(i1,i2))) ;

		xv[3] = dum*cos(theta) ;
    xv[4] = dum*sin(theta) ;
		theta = 2.0*pi*SimRan1(i1,i2);

    // dum = vexp*sqrt(-log(SimRan1(i1,i2))) ;
    dum = sigma*sqrt(-2.0*log(SimRan1(i1,i2))) ;

		xv[5] = dum*cos(theta) ;
		dum = 0.5*mprot*(xv[3]*xv[3] + xv[4]*xv[4] + xv[5]*xv[5]); // KE of Hbar
	}
}



//--------- random number generator, returns between 0 and 1
double SimRan1(long int& i1, long int& i2) {
	long int k, j1;
	i1 = 16807*(i1%127773) - 2836*(i1/127773); if(i1 <= 0) i1+= 2147483647;
	i2 = 40014*(i2%53668) - 12211*(i2/53668);  if(i2 <= 0) i2+= 2147483563;
	k = i2/ndiv2; j1 = iv1[k]; iv1[k] = i1;
	return ((j1 - 0.5)*4.656612875245797E-10);
}

//--------- Calculate the energy of the Hbar, Kinetic and (Mag) Potential
double Energy(double xv[]){

	double T, U, bft[3], bmag;

  // kinetic term
  T = 0.5*mprot*(xv[3]*xv[3]+xv[4]*xv[4]+xv[5]*xv[5]);

  // potential term
	bft[0] = 0.0; bft[1] = 0.0; bft[2] = 0.0;
	Bfield(xv,bft);
	bmag = sqrt(bft[0]*bft[0]+bft[1]*bft[1]+bft[2]*bft[2]) ;
	U = magmom0*(bmag-bfmin);  // not a dot product?
	return T + U ;
}


//--------- this computes the magnetic field for the ALPHA trap
//          includes the mirror coils, the uniform field, and the octupole field
void Bfield(double r[], double bf[]) {
  if (strcmp(potential_option, "ho") == 0) { // Harmonic Oscillator
    bf[0] = 1.54*(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]) / (radtrp2);
    bf[1] = 0;
    bf[2] = 0;
  } else {
  	double bf_temp[3] ;
  	// MIRROR COILS
    // use this for simplified calculation of solenoid:
  	Bsolesimp(r,bf_temp);
    // use this for more accurate calculation of solenoid:
  	// Bsoles(adi,zpos0,zposf,bcent,r,bfp);
  	for (int j = 0 ; j < 3 ; j++) {
      bf[j] = bf_temp[j];
    }
    if (strcmp(potential_option, "oct") == 0) { // Octupole
  	  Boct(r,bf_temp);
    } else if (strcmp(potential_option, "quad") == 0) { // Quadrupole
      Bquad(r, bf_temp);
    }
    for (int j = 0 ; j < 3 ; j++) {
      bf[j] += bf_temp[j];
    }

  	// UNIFORM FIELD
  	bf[2] += bunif;
  }
}


//--------- Magnetic field from solenoid
/* It approximates the solenoid as 2 loops separated by the distance
   given by 2 point Gauss-Legendre integration +-L/sqrt(3)
   It uses the approximation 5.39 of Jackson for the vector potential
   to compute the magnetic field B = curl(A) */
void Bsolesimp(double r[], double bf[]) {
  double bfp[3], bfold[3], r2a, s2, s, r2b, s2r, r2rs2r, r3a;

  s2 = r[0]*r[0] + r[1]*r[1] ; s = sqrt(s2) ;

  if(s > (1.e-7*radis)) {
    //compute B from mirr1, loop 1 using ansatz
    r2a = sqrt(s2 + (r[2]-zpos11)*(r[2]-zpos11) + radis2 - 2*lam1*s*radis) ;
    r2b = sqrt(s2 + (r[2]-zpos11)*(r[2]-zpos11) + radis2 + 2*lam1*s*radis) ;
    bfp[1] = imu0a1*(r[2]-zpos11)*(1/(r2a*r2a*r2a) - 1/(r2b*r2b*r2b))/(s*radis*2*lam1) ;
    bfp[2] = imu0a1*(1/r2a - 1/r2b - s*(s-lam1*radis)/(r2a*r2a*r2a) + s*(s+lam1*radis)/(r2b*r2b*r2b))/(s*radis*2*lam1) ;

    //compute B from mirr1, loop 2 using ansatz
    r2a = sqrt(s2 + (r[2]-zpos12)*(r[2]-zpos12) + radis2 - 2*lam1*s*radis) ;
    r2b = sqrt(s2 + (r[2]-zpos12)*(r[2]-zpos12) + radis2 + 2*lam1*s*radis) ;
    bfold[1] = imu0a1*(r[2]-zpos12)*(1/(r2a*r2a*r2a) - 1/(r2b*r2b*r2b))/(s*radis*2*lam1) ;
    bfold[2] = imu0a1*(1/r2a - 1/r2b - s*(s-lam1*radis)/(r2a*r2a*r2a) + s*(s+lam1*radis)/(r2b*r2b*r2b))/(s*radis*2*lam1) ;

    //compute the B for the solenoid 1
    for(int j = 1 ; j < 3 ; j++) {
      bf[j] = 0.5*(bfold[j]+bfp[j]);
    }

    //compute B from mirr2, loop 1 using ansatz
    r2a = sqrt(s2 + (r[2]-zpos21)*(r[2]-zpos21) + radis2 - 2*lam2*s*radis) ;
    r2b = sqrt(s2 + (r[2]-zpos21)*(r[2]-zpos21) + radis2 + 2*lam2*s*radis) ;
    bfp[1] = imu0a2*(r[2]-zpos21)*(1/(r2a*r2a*r2a) - 1/(r2b*r2b*r2b))/(s*radis*2*lam2) ;
    bfp[2] = imu0a2*(1/r2a - 1/r2b - s*(s-lam2*radis)/(r2a*r2a*r2a) + s*(s+lam2*radis)/(r2b*r2b*r2b))/(s*radis*2*lam2) ;

    //compute B from mirr2, loop 2 using ansatz
    r2a = sqrt(s2 + (r[2]-zpos22)*(r[2]-zpos22) + radis2 - 2*lam2*s*radis) ;
    r2b = sqrt(s2 + (r[2]-zpos22)*(r[2]-zpos22) + radis2 + 2*lam2*s*radis) ;
    bfold[1] = imu0a2*(r[2]-zpos22)*(1/(r2a*r2a*r2a) - 1/(r2b*r2b*r2b))/(s*radis*2*lam2) ;
    bfold[2] = imu0a2*(1/r2a - 1/r2b - s*(s-lam2*radis)/(r2a*r2a*r2a) + s*(s+lam2*radis)/(r2b*r2b*r2b))/(s*radis*2*lam2) ;

    //compute the B for the solenoid 2
    for(int j = 1 ; j < 3 ; j++) {
      bf[j] += 0.5*(bfold[j]+bfp[j]);
    }

    //convert B_s/s to B_x and B_y
    bf[0] = bf[1]*r[0] ; bf[1] *= r[1] ;
  } else {
    //compute B from mirr1, loop 1 using ansatz
    r2a = 1/(s2 + (r[2]-zpos11)*(r[2]-zpos11) + radis2) ;
    s2r = s2*r2a ; r2rs2r = radis2*r2a*s2r ;
    r3a = imu0a1*r2a*sqrt(r2a) ;
    bfp[1] = (r[2]-zpos11)*(3+17.5*lam1*lam1*r2rs2r)*r2a*r3a ;
    bfp[2] = (2-3*s2r+10*lam1*lam1*r2rs2r*(1-1.75*s2r))*r3a ;


    //compute B from mirr1, loop 2 using ansatz
    r2a = 1/(s2 + (r[2]-zpos12)*(r[2]-zpos12) + radis2) ;
    s2r = s2*r2a ; r2rs2r = radis2*r2a*s2r ;
    r3a = imu0a1*r2a*sqrt(r2a) ;
    bfold[1] = (r[2]-zpos12)*(3+17.5*lam1*lam1*r2rs2r)*r2a*r3a ;
    bfold[2] = (2-3*s2r+10*lam1*lam1*r2rs2r*(1-1.75*s2r))*r3a ;

    //compute the B for the solenoid 1
    for (int j = 1 ; j < 3 ; j++) {
      bf[j] = 0.5*(bfold[j]+bfp[j]);
    }


    //compute B from mirr2, loop 1 using ansatz
    r2a = 1/(s2 + (r[2]-zpos21)*(r[2]-zpos21) + radis2) ;
    s2r = s2*r2a ; r2rs2r = radis2*r2a*s2r ;
    r3a = imu0a2*r2a*sqrt(r2a) ;
    bfp[1] = (r[2]-zpos21)*(3+17.5*lam1*lam1*r2rs2r)*r2a*r3a ;
    bfp[2] = (2-3*s2r+10*lam1*lam1*r2rs2r*(1-1.75*s2r))*r3a ;

    //compute B from mirr2, loop 2 using ansatz
    r2a = 1/(s2 + (r[2]-zpos22)*(r[2]-zpos22) + radis2) ;
    s2r = s2*r2a ; r2rs2r = radis2*r2a*s2r ;
    r3a = imu0a2*r2a*sqrt(r2a) ;
    bfold[1] = (r[2]-zpos22)*(3+17.5*lam1*lam1*r2rs2r)*r2a*r3a ;
    bfold[2] = (2-3*s2r+10*lam1*lam1*r2rs2r*(1-1.75*s2r))*r3a ;

    //compute the B for the solenoid 2
    for(int j = 1 ; j < 3 ; j++) {
      bf[j] += 0.5*(bfold[j]+bfp[j]);
    }

    //convert B_s/s to B_x and B_y
    bf[0] = bf[1]*r[0] ; bf[1] *= r[1] ;
  }
}

//--------- Magnetic field from OCTOPOLE
/* it uses a fit to the functional form for the vector potential
   where the A_z step is treated as the difference of 2 error functions
   for z > zmid, the largest error in |B_dirk-B_fit| ~ 1%
   ... AND ...
   it does not include the effect from the leads
   imu1 is the I mu_0/(4 pi) */
void Boct(double r[], double bf[]) {
  double x1, x2, d0, d1, d2, d3, d4, f40, f41, f42, f43, f44, f45, f46, s2, x4,
         x2y2, y4, dum;

  x2 = (r[2]-zmid-z1)/dz1 ; x1 = (r[2]-zmid+z1)/dz1 ;

  if((x2 > 4.5) || (x1 < -4.5)) {
    bf[0] = 0.0 ; bf[1] = 0.0 ; bf[2] = 0.0 ;
  } else {
    if ((x2 < -4.5) && (x1 > 4.5)) {
      bf[0] = r[1]*4*2*imu1*(r[1]*r[1]-3*r[0]*r[0]) ;
      bf[1] =-r[0]*4*2*imu1*(r[0]*r[0]-3*r[1]*r[1]) ;
      bf[2] = 0.0 ;
    } else {
      if (x1 <= 4.5) {
        f40 = imu1*(2.0 - Erfcc(x1)) ;
        dum = exp(-x1*x1)*imu1 ;
        f41 = dum*co1 ;
        f42 =-dum*x1*co2 ;
        f43 = dum*(x1*x1-0.5)*co3 ;
        f44 = dum*x1*(1.5-x1*x1)*co4 ;
        f45 = dum*(.75-3.0*x1*x1+x1*x1*x1*x1)*co5 ;
        f46 =-dum*x1*(3.75-5.0*x1*x1+x1*x1*x1*x1)*co6 ;

        s2 = r[0]*r[0] + r[1]*r[1] ;
        d0 = f40 - 0.05*f42*s2 + f44*s2*s2/960;
        d1 =-0.1*(f42 - f44*s2/24 + f46*s2*s2/1344) ;
        d2 =-0.1*f42 + f44*s2/240 ;
        d3 =-0.1*(f41 - f43*s2/24 + f45*s2*s2/1344) ;
        d4 = -(-f43/120 + f45*s2/3360);

        x4 = r[0]*r[0]*r[0]*r[0] ;
        x2y2 = r[0]*r[0]*r[1]*r[1] ;
        y4 = r[1]*r[1]*r[1]*r[1] ;

        bf[0] = r[1]*(4*d0*(r[1]*r[1]-3*r[0]*r[0]) + d2*(x4-6*x2y2+y4)
             -d1*(y4-10*x2y2+5*x4)) ;
        bf[1] =-r[0]*(4*d0*(r[0]*r[0]-3*r[1]*r[1]) + d2*(y4-6*x2y2+x4)
             -d1*(x4-10*x2y2+5*y4)) ;
        bf[2] = 4*(10*d3 + d4*s2)*r[0]*r[1]*(r[0]*r[0]-r[1]*r[1]) ;
      } else {
        f40 = imu1*Erfcc(x2) ;
        dum = imu1*exp(-x2*x2) ;
        f41 =-dum*co1 ;
        f42 = dum*x2*co2 ;
        f43 =-dum*(x2*x2-0.5)*co3 ;
        f44 =-dum*x2*(1.5-x2*x2)*co4 ;
        f45 =-dum*(0.75-3.0*x2*x2+x2*x2*x2*x2)*co5 ;
        f46 = dum*x2*(3.75-5.0*x2*x2+x2*x2*x2*x2)*co6 ;

        s2 = r[0]*r[0] + r[1]*r[1] ;
        d0 = f40 - 0.05*f42*s2 + f44*s2*s2/960;
        d1 =-0.1*(f42 - f44*s2/24 + f46*s2*s2/1344) ;
        d2 =-0.1*f42 + f44*s2/240 ;
        d3 =-0.1*(f41 - f43*s2/24 + f45*s2*s2/1344) ;
        d4 = -(-f43/120 + f45*s2/3360);

        x4 = r[0]*r[0]*r[0]*r[0] ;
        x2y2 = r[0]*r[0]*r[1]*r[1] ;
        y4 = r[1]*r[1]*r[1]*r[1] ;

        bf[0] = r[1]*(4*d0*(r[1]*r[1]-3*r[0]*r[0]) + d2*(x4-6*x2y2+y4)
             -d1*(y4-10*x2y2+5*x4)) ;
        bf[1] =-r[0]*(4*d0*(r[0]*r[0]-3*r[1]*r[1]) + d2*(y4-6*x2y2+x4)
             -d1*(x4-10*x2y2+5*y4)) ;
        bf[2] = 4*(10*d3 + d4*s2)*r[0]*r[1]*(r[0]*r[0]-r[1]*r[1]) ;
      }
    }
  }
}


void Bquad(double r[], double bf[]) {
  bf[0] = imu2*r[1];
  bf[1] = imu2*r[0];
  bf[2] = 0.0 ;
}

void Boct_edgeless(double r[], double bf[]) {
    bf[0] = r[1]*4*2*imu1*(r[1]*r[1]-3*r[0]*r[0]) ;
    bf[1] =-r[0]*4*2*imu1*(r[0]*r[0]-3*r[1]*r[1]) ;
    bf[2] = 0.0 ;
}



//--------- Function that calculates the derivative of the 1st order ODE
void Dydt(double t, double xv[], double dxvdt[]){
	double frcr, frcz, bf[3], frc[3];

  for (int j = 0 ; j < 3 ; j++) {
		dxvdt[j] = xv[j+3];
	}

	// Calculate magnetic force with updated current, etc.
	ForceB(xv,frc);

	// Put in effect of trapping fields.
  for (int j = 0 ; j < 3 ; j++) {
		dxvdt[j+3] = frc[j] / mprot;
	}
}


// ------- loop continue condition in main()
bool ContinueLoop(double tim, double tfin, double xv[3]) {
  // original script
  /*if(tim < tfin){
    if ((xv[2] > 0.0) && (xv[2] < 0.04516)){
      if ((xv[0]*xv[0]+xv[1]*xv[1]) < 0.016800 * 0.016800){ l_cond = 1; }
    }
    if ((xv[2] >= 0.04516) && (xv[2] <= 0.32036)){
      if ((xv[0]*xv[0]+xv[1]*xv[1]) < 0.022275 * 0.022275){ l_cond = 1; }
    }
    if ((xv[2] > 0.32036) && (xv[2] < 0.37848)){
      if ((xv[0]*xv[0]+xv[1]*xv[1]) < 0.016800){ l_cond = 1; }
    }
  }
  zmid = 0.183385;

  radtrp = 0.022275;   // Trap radius (maximum electrode radius)
  sm_rad = 0.016800;   // Trap radius (small electrode radius) */

  if (tim < tfin) {
    if ((xv[2] > end1) && (xv[2] < end5)) { // within ends of trap
      if ((xv[2] < end2) || (xv[2] > end4)) { // small trap
        if ((xv[0]*xv[0]+xv[1]*xv[1]) < sm_rad2)
          return true;
      } else if (xv[2] < end3) { // mid trap
        if ((xv[0]*xv[0]+xv[1]*xv[1]) < md_rad2)
          return true;
      } else {
        if ((xv[0]*xv[0]+xv[1]*xv[1]) < radtrp2)
          return true;
      }
    }
  }
  return false;
}

// ------- Symplectic Stepper!!
void Symplec(double tim, double dtim, double xv[]){
	double dum1, dum2, frcr, frcz, bf[3], vec[3], frc[3];

	// x(t - dt/2) to x(t + dt/2) using v(t)
	for (int j = 0 ; j < 3 ; j++) {
		xv[j] += (dtim*xv[j+3]);
	}

	//Calculate magnetic force

  // v(t) to v(t + dt) using x(t + dt/2)
	ForceB(xv,frc);
  for (int j = 0 ; j < 3 ; j++) {
		xv[j+3] += dtim*(frc[j] / mprot);
	}
}


// -------- This routine computes the MAGNETIC force felt by the Hbar.
void ForceB(double xv[], double frc[]) {
	double dx, bp, bm, rp[3], rm[3], bf[3];

	//dx is needed for numerical derivative of |B| to compute force on Hbar
	dx = radtrp*1.e-5 ;

  for (int i = 0 ; i < 3 ; i++) {
    rp[i] = xv[i];
    rm[i] = xv[i];
  }

	for (int i = 0 ; i < 3 ; i++) {
		rp[i] += dx/2 ; rm[i] -= dx/2 ;

		Bfield(rp,bf);
		bp = sqrt(Dot(bf, bf));
		Bfield(rm,bf);
		bm = sqrt(Dot(bf, bf));

		frc[i] = -magmom0*(bp - bm)/dx;

		rp[i] = xv[i]; rm[i] = xv[i];
	}
}

std::string DateString() {
  char buff[20];
  time_t now = time(NULL);
  strftime(buff, 20, "%Y%m%d_%H%M%S", localtime(&now));
  return (std::string) buff;
}

// ------- Writes metaparameters –– the parameters shared across all trajectories
void WriteMetaParams(std::string name, int numtraj, double tfin, double dtim, int seed) {
  std::ofstream metaparams;
  metaparams.open(("metaparams_" + name + ".csv").c_str(), std::ios_base::app);
  metaparams << "seed,numtraj,tfin,dtim\n";
  metaparams << seed << "," << numtraj << "," << tfin << "," << dtim << "\n";
  metaparams.close();
}


void WriteCrossings(std::string name, std::vector<double> cross_t, std::vector<double> cross_Tz) { // std::vector<double> cross_x, std::vector<double> cross_y, std::vector<double> cross_vx, std::vector<double> cross_vy) {
  std::ofstream crossings ((name + "_crossings.csv").c_str());
  crossings.is_open();
  crossings << "cross_t,cross_dt,cross_Tz\n";
  crossings << cross_t[0] << "," << 0 << ',' << cross_Tz[0] << "," << 0 << "\n";
  for (int i = 1; i < cross_t.size(); i++) {
    crossings << cross_t[0] << "," << (cross_t[i] - cross_t[i-1]) << "," << cross_Tz[i] << "\n";
  }
  crossings.close();
}

long double average(std::vector<double> cross_Tz) {
  long double avg;
  for (int i = 0; i < cross_Tz.size(); i++) {
    avg += cross_Tz[i];
  }
  avg /= cross_Tz.size();
  return avg;
}

long double variance(std::vector<double> cross_Tz, long double average) {
  long double variance;
  for (int i = 0; i < cross_Tz.size(); i++) {
    variance += (cross_Tz[i] - average)*(cross_Tz[i] - average);
  }
  variance /= cross_Tz.size();
  return variance;
}
