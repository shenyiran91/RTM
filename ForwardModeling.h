#ifndef FORWARDMODELING_H
#define FORWARDMODELING_H

#include <array>
#include <vector>
#include <random>
#include "Source.h"

const double R = 0.0000000001; // relative reflection
const double pi{3.1415926};

const double m_mean = 0.0; // data white noise
const double m_stddev = 0.00;
//const unsigned int seed = time(0);
class ForwardModeling
{
private:
    int m_nx, m_nz;   // size of the computation area
    double m_dx, m_dz; // spatial interval of the computation area
    int m_nt;         // number of steps of forward modeling
    double m_dt;       // forward modeling time interval
    int m_nb;         // width of the pml boundary area
    bool m_back_fd, save;
    
    std::array<double, 5> coeff{-205. / 72, 8. / 5, -1. / 5, 8. / 315, -1. / 560};
    std::array<double, 5> coeff_1fd{ 0.0 , -4.0 / 5.0, 1.0 / 5.0, -4.0 / 105.0, 1.0 / 280.0}; 
    //std::array<double, 4> coeff{-49. / 18, 3. / 2, -3. / 20, 1. / 90,};
    //std::array<double, 3> coeff{-5./ 2, 4./ 3, -1./12};
    //std::array<double, 2> coeff{-2., 1.};
    int fd_order, half_fd_order, fd_order_half;  // finite difference order (half)
    int m_global_nx, m_global_nz; // nx + 2 * nb + 2 * 4 (fd order)
    int start_x, start_z;
    
    Source m_source; // source wavelet

    std::vector<double> p_old, p_cur, p_new, p_bc; // wavefield
    
    int n_depth, n_depth_1, n_depth_2;
    std::vector<double> v;          // velocity field
    std::vector<double> v2;          // velocity field
    std::vector<double> v3;          // velocity field    
    std::vector<double> a;          // inverse density field
    std::vector<double> a2;
    std::vector<double> a3;
    //std::vector<double> a4;
    std::vector<double> b;          // bulk modulus
    std::vector<double> b2;
    std::vector<double> b3;
    //std::vector<double> b4;
    std::vector<double> force, boundary_force;      // force term
    std::vector<double> zeta_x, zeta_z;     // damping profile
    std::vector<double> phi_cur_x, phi_cur_z; // auxiliary variables
    std::vector<double> phi_new_x, phi_new_z; // auxiliary variables
    std::vector<double> bcp_new;     // boundary data
    std::vector<double> bcstep;     // boundary data collect


    double m_dx2, m_dz2, m_dt2;

    int getIdx(int ix, int iy);

public:
    ForwardModeling(int nx, int nz, double dx, double dz, int nt, double dt,
                    int nb, bool back_fd, const Source &source, const std::vector <double> &parameters);
    double g_random_nosie(double mean, double std);
    void setForce(const double &t);
    void setbackForce(const double &t);
    void setDampingProfile();
    double laplacian(const int &ix, const int &iz);
    double central_FD_first_derivative(const int &ix, const int &iz);
    double laplacian2(const int &ix, const int &iz);
    void updatePhi(const int &ix, const int &iz);
    bool isBoundary(const int &ix, const int &iz);
    bool isTopBoundary(const int &ix, const int &iz);
    std::vector<double> step(double t);
    void runstep();
};

#endif
