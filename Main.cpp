#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <random>
#include "ForwardModeling.h"
#include "Source.h"
#include <math.h>
#include <fstream>
#include <vector>
#include <stdexcept>
//#include <limits>


// Defind some parameters for test
// Later, implement reading parameters from file
const int nx = 100;    //depth
const int nz = 100; //200 length
const double dx = 10;
const double dz = 10;
const double nt = 1000;
const double dt = 0.001;
const int nb = 50; // boundary width
const double source_a = 5.0;
const double source_f0 = 30;
const double source_t0 = 0.02;
const int source_locx = 18; //18
const int source_locz = 50; //100

int main()
 {
  std::vector<double> parameters;
  parameters.push_back(4.0e+06);  // v^2 4000000 3.99917e+06
  parameters.push_back(0.5);    // 0.5   0.562056
  parameters.push_back(600);    // 600.    760.053

  Source source = Source(source_locx, source_locz, source_a, source_f0, source_t0);

  //ForwardModeling forward = ForwardModeling(nx, nz, dx, dz, nt, dt, nb, false, source, parameters);
  //forward.runstep();

  ForwardModeling backward = ForwardModeling(nx, nz, dx, dz, nt, dt, nb, true, source, parameters);
  backward.runstep();

 // std::vector<double> bcp_new(100,0);
 // std::vector<double> force( (100 + 2 * 50 + 8) * (100 + 2 * 50 + 8), 0);
 // for (size_t ii = 20; ii < 21; ++ii) {
 //     double t  = ii*10;
 // 
 //   // Open the file in binary mode
 //   std::ifstream infile("Boundaryfield.dat", std::ios::binary);
 //   if (!infile) {
 //     throw std::runtime_error("Cannot open Boundaryfield.dat");
 //   }
 //   // Calculate the position to seek to based on time t
 //   std::streampos vector_position = t * bcp_new.size() * sizeof(double);

 //   // Seek to the position
 //   infile.seekg(vector_position);
 //   if (!infile) {
 //     throw std::runtime_error("Error seeking to position in Boundaryfield.dat");
 //   }

 //   // Read the data into the force vector
 //   infile.read(reinterpret_cast<char*>(bcp_new.data()), bcp_new.size() * sizeof(double));
 //   if (!infile) {
 //     throw std::runtime_error("Error reading from Boundaryfield.dat");
 //   }

 //   // Close the file
 //   infile.close();
 //     
 //   std::fill(force.begin(), force.end(), 0);
 //   for (size_t i = 0; i < bcp_new.size(); ++i) {
 //     force[i+54] = bcp_new[i];
 //     printf("force[%d] = %f\n", i+54, force[i+54]);
 //   }
 // }
 //   // Create a .dat file with the vector data
 //   std::ofstream dataFile("data.dat");
 //   for (size_t i = 0; i < bcp_new.size(); ++i) {
 //       dataFile << i << " " << bcp_new[i] << std::endl;  // Write index and value
 //   }
 //   dataFile.close();

 //   // Create a GNUplot script
 //   std::ofstream plotScript("plot.plt");
 //   plotScript << "set title 'Vector Plot'\n";
 //   plotScript << "set xlabel 'Index'\n";
 //   plotScript << "set ylabel 'Value'\n";
 //   plotScript << "plot 'data.dat' using 1:2 with linespoints title 'Data'\n";
 //   plotScript.close();

 //   // Run the plot script with GNUplot
 //   system("gnuplot -p plot.plt");


 // parameters.clear();
  return 0;
}
