#include "ForwardModeling.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>
#include <cassert>
#include <ctime>
#include <iterator>
#include <fstream>
#include <vector>
#include <stdexcept>

ForwardModeling::ForwardModeling(int nx, int nz, double dx, double dz, int nt,
                                 double dt, int nb, bool back_fd, const Source &source, const std::vector <double> &parameters)
    : m_nx{nx}, m_nz{nz}, m_dx{dx}, m_dz{dz}, m_nt{nt}, m_dt{dt}, m_nb{nb}, m_back_fd{back_fd},
      fd_order{(static_cast<int>(coeff.size()) - 1) * 2},
      fd_order_half{(static_cast<int>(coeff.size()) - 1)},
      m_global_nx{nx + 2 * nb + fd_order},
      m_global_nz{nz + 2 * nb + fd_order},
      m_source{source},
      p_old(m_global_nx * m_global_nz, 0),
      p_cur(m_global_nx * m_global_nz, 0),
      p_new(m_global_nx * m_global_nz, 0),
      p_bc(m_global_nx, 0),
      n_depth{(int)floor(parameters[2]/m_dx)}, 
      v(m_global_nx * m_global_nz, 2000.),
      a(m_global_nx * ((n_depth) + nb + fd_order_half), 1000.),
      a2(m_global_nx * ((nz-n_depth) + nb + fd_order_half), 1000. * (1 + parameters[1])),
      b(m_global_nx * ((n_depth) + nb + fd_order_half), parameters[0]/ 1000.),
      b2(m_global_nx * ((nz-n_depth) + nb + fd_order_half), parameters[0]/ 1000. / (1 + parameters[1])),
      
      force(m_global_nx * m_global_nz, 0),
      boundary_force(m_nz, 0),
      zeta_x(m_global_nx * m_global_nz, 0),
      zeta_z(m_global_nx * m_global_nz, 0),
      phi_cur_x(m_global_nx * m_global_nz, 0),
      phi_cur_z(m_global_nx * m_global_nz, 0),
      phi_new_x(m_global_nx * m_global_nz, 0),
      phi_new_z(m_global_nx * m_global_nz, 0),
      bcp_new(m_nz,0),
      bcstep(m_nz,0)
{
    half_fd_order = fd_order / 2;
    start_x = half_fd_order + m_nb;
    start_z = half_fd_order + m_nb;
    m_dx2 = m_dx * m_dx;
    m_dz2 = m_dz * m_dz;
    m_dt2 = m_dt * m_dt;

    save = true;


    a.insert(a.end(), a2.begin(), a2.end());
    b.insert(b.end(), b2.begin(), b2.end());
   
    setDampingProfile();
}

/*************************************************
 * Get the index in 1D vector given 2D location
 * @param ix the horizontal index of the point.
 * @param iz the vertical index of the point.
 * @return the index in 1D vector.
 *************************************************/
int ForwardModeling::getIdx(int ix, int iz)
{
    assert(ix >= 0);
    assert(ix < m_global_nx);
    assert(iz >= 0);
    assert(iz < m_global_nz);

    return ix * m_global_nz + iz;
}

/*************************************************
 * Set the force term for the current time using
 * ricker wavelet
 * @param the current time
 *************************************************/
void ForwardModeling::setForce(const double &t)
{
    std::fill(force.begin(), force.end(), 0);
    force[getIdx(m_source.getX() + start_x, m_source.getZ() + start_z)] =
        m_source.rickerWavelet(t);
}

/*************************************************
 * Set the force term for the current time using
 * gather
 * @param the current time
 *************************************************/
void ForwardModeling::setbackForce(const double &t)
{
    // Open the file in binary mode
    std::ifstream infile("Boundaryfield.dat", std::ios::binary);
    if (!infile) {
        throw std::runtime_error("Cannot open Boundaryfield.dat");
    }
    // Calculate the position to seek to based on time t
    std::streampos vector_position = (m_nt - (t * 1/m_dt) - 1) * bcp_new.size() * sizeof(double);

    // Seek to the position
    infile.seekg(vector_position);
    if (!infile) {
        throw std::runtime_error("Error seeking to position in Boundaryfield.dat");
    }

    // Read the data into the force vector
    infile.read(reinterpret_cast<char*>(boundary_force.data()), boundary_force.size() * sizeof(double));
    if (!infile) {
        throw std::runtime_error("Error reading from Boundaryfield.dat");
    }

    if(t >= 0.65){
        std::fill(boundary_force.begin(), boundary_force.end(), 0);
    }

    std::fill(force.begin(), force.end(), 0);
    for (size_t i = 0; i < boundary_force.size(); ++i) {
        force[ ( ( start_x + 1 ) * m_global_nz ) + i+ start_z + 1] = boundary_force[i];
        //if( t == 0.799){
        //    printf("force[%d] = %f\n", ( (start_x + 1) * m_global_nz) + i + start_z + 1, force[( (start_x+1) * m_global_nz) + i+ start_z + 1]);
        //}
    }

    //std::ofstream BCf{"force.dat", std::ios::binary | std::ios::app};
    //BCf.write((char *)&force[0],force.size() * sizeof(double));

    //std::ofstream BCft{"force_test_longer.dat", std::ios::binary | std::ios::app};
    //BCft.write((char *)&boundary_force[0],boundary_force.size() * sizeof(double));
    
    // Close the file
    infile.close();
}

/************************************************
 * Set the damping profile for each point
 * zeta = \frac{c}{L}log(\frac{1}{R}) * 
 *        (\frac{|x-a|}{L} - 
 *         \frac{sin(\frac{2\pi|x-a|}{L})}{2\pi})
 ************************************************/
void ForwardModeling::setDampingProfile()
{
    // left boundary
    for (int ix{0}; ix < start_x; ix++)
        for (int iz{0}; iz < m_global_nz; iz++)
        {
            double idx = getIdx(ix, iz);
            double zeta = v[idx] / (start_x * m_dx) * (-log(R));
            //double zeta = sqrt (a[idx] * b[idx]) / (start_x * m_dx) * (-log(R));
            double r_length = 1. * (start_x - ix) / start_x;
            zeta_x[idx] = zeta * (r_length - sin(2. * pi * r_length) / 2. / pi);
        }
    // right boundary
    for (int ix{m_global_nx - start_x}; ix < m_global_nx; ix++)
        for (int iz{0}; iz < m_global_nz; iz++)
        {
            double idx = getIdx(ix, iz);
            double zeta = v[idx] / (start_x * m_dx) * (-log(R));
            //double zeta = sqrt (a[idx] * b[idx]) / (start_x * m_dx) * (-log(R));
            double r_length = 1. * (ix + start_x - m_global_nx) / start_x;
            zeta_x[idx] = zeta * (r_length - sin(2. * pi * r_length) / 2. / pi);
        }
    // top boundary
    for (int ix{0}; ix < m_global_nx; ix++)
        for (int iz{0}; iz < start_z; iz++)
        {
            double idx = getIdx(ix, iz);
            double zeta = v[idx] / (start_z * m_dz) * (-log(R));
            //double zeta = sqrt (a[idx] * b[idx]) / (start_z * m_dz) * (-log(R));
            double r_length = 1. * (start_z - iz) / start_z;
            zeta_z[idx] = zeta * (r_length - sin(2. * pi * r_length) / 2. / pi);
        }
    // bottom boundary
    for (int ix{0}; ix < m_global_nx; ix++)
        for (int iz{m_global_nz - start_z}; iz < m_global_nz; iz++)
        {
            double idx = getIdx(ix, iz);
            double zeta = v[idx] / (start_z * m_dz) * (-log(R));
            //double zeta = sqrt (a[idx] * b[idx])  / (start_z * m_dz) * (-log(R));
            double r_length = 1. * (iz + start_z - m_global_nz) / start_z;
            zeta_z[idx] = zeta * (r_length - sin(2. * pi * r_length) / 2. / pi);
        }
}

bool ForwardModeling::isBoundary(const int &ix, const int &iz)
{
    if (ix > start_x && ix < m_nx + start_x && iz > start_z && iz < m_nz + start_z)
        return false;
    else
        return true;
}


bool ForwardModeling::isTopBoundary(const int &ix, const int &iz)
{
    if (ix == start_x+1 && iz > start_z && iz < m_nz + start_z )
    //if (ix == start_x+1)
        return true;
    else
        return false;
}


/*************************************************
 * Laplacian operator
 * \Delta f = \sum_{i=-4}^4 c[abs(i)]*p[ix-i][iz] / dx2
 *          + \sum_{i=-4}^4 c[abs(i)]*p[ix][iz-i] / dz2
 *************************************************/
double ForwardModeling::laplacian(const int &ix, const int &iz)
{
    auto fd_sum = p_cur[getIdx(ix, iz)] * coeff[0] / m_dx2 +
                  p_cur[getIdx(ix, iz)] * coeff[0] / m_dz2;
    for (size_t i{1}; i < coeff.size(); i++)
    {
        fd_sum += (p_cur[getIdx(ix - i, iz)] +
                   p_cur[getIdx(ix + i, iz)]) *
                      coeff[i] / m_dx2 +
                  (p_cur[getIdx(ix, iz - i)] +
                   p_cur[getIdx(ix, iz + i)]) *
                      coeff[i] / m_dz2;
    }
    return fd_sum;

}
/*************************************************
 * Central finite difference first derivative
 * \partial f / \partial x = coeff_1st * (f(x+i*dx) - f(x-i*dx)) / dx
 *************************************************/
double ForwardModeling::central_FD_first_derivative(const int &ix, const int &iz)
{
    double db_x = 0.0; //(b[getIdx(ix + 1, iz)] - b[getIdx(ix - 1, iz)]) / 2;
    double dp_x = 0.0; //(p_cur[getIdx(ix + 1, iz)] - p_cur[getIdx(ix - 1, iz)]) / 2;
    double db_y = 0.0; //(b[getIdx(ix, iz + 1)] - b[getIdx(ix , iz - 1)]) / 2;
    double dp_y = 0.0; //(p_cur[getIdx(ix , iz + 1)] - p_cur[getIdx(ix, iz - 1)]) / 2;;
    for (size_t i{1}; i < coeff_1fd.size(); i++)
    {
        db_x += (b[getIdx(ix  + i, iz)]  - b[getIdx(ix  - i, iz)]) * coeff_1fd[i];
        dp_x += (p_cur[getIdx(ix + i, iz)] - p_cur[getIdx(ix - i, iz)]) * coeff_1fd[i];
        db_y += (b[getIdx(ix, iz + i)] - b[getIdx(ix , iz - i)]) * coeff_1fd[i];
        dp_y += (p_cur[getIdx(ix, iz + i)] - p_cur[getIdx(ix, iz - i)])* coeff_1fd[i];
    }
    double fd = (db_x  * dp_x)/ m_dx  + (db_y * dp_y)/ m_dz;
    return fd;
}


///*************************************************
// * Laplacian operator 2  
// * (d/dx b * (d p_cur)/dx) = (db(i)/dx * dp/dx + db(i)/dz * dp/dz) + b(i) * (d^2 p/dx^2 + d^2 p/dz^2)
// *        \Delta f         = \first term + b(i) * laplacian
// *                         = 1st_derivative + b(i) * laplacian
///*************************************************/
//double ForwardModeling::laplacian2(const int &ix, const int &iz)
//{
//    double first_term = central_FD_first_derivative(ix, iz);
//    double fd_sum = first_term + b[getIdx(ix, iz)] * laplacian(ix, iz);
//    return fd_sum;
//}



/*************************************************
 * Laplacian operator 2
 * \Delta f = \sum_{i=-4}^4 c[abs(i)]*p[ix-i][iz] *(avg_b(i))/ dx2
 *          + \sum_{i=-4}^4 c[abs(i)]*p[ix][iz-i] *(avg_b(i)) / dz2
 *************************************************/
double ForwardModeling::laplacian2(const int &ix, const int &iz)
{
    auto fd_x1 = 0.5 * p_cur[getIdx(ix, iz)] * coeff[0] * (b[getIdx(ix + 1, iz)] + b[getIdx(ix , iz)]) /2;
    auto fd_x2 = 0.5 * p_cur[getIdx(ix, iz)] * coeff[0] * (b[getIdx(ix , iz)] + b[getIdx(ix - 1, iz)]) /2;
    auto fd_z1 = 0.5 * p_cur[getIdx(ix, iz)] * coeff[0] * (b[getIdx(ix , iz + 1)] + b[getIdx(ix , iz)]) /2;
    auto fd_z2 = 0.5 * p_cur[getIdx(ix, iz)] * coeff[0] * (b[getIdx(ix , iz)] + b[getIdx(ix, iz - 1)]) /2;
    for (size_t i{1}; i < coeff.size(); i++)
    {
        fd_x1 += p_cur[getIdx(ix+i, iz)] * coeff[i] * (b[getIdx(ix + i, iz)] + b[getIdx(ix , iz)]) /2;
        fd_x2 += p_cur[getIdx(ix-i, iz)] * coeff[i] * (b[getIdx(ix , iz)] + b[getIdx(ix - i, iz)]) /2;
        fd_z1 += p_cur[getIdx(ix, iz+i)] * coeff[i] * (b[getIdx(ix, iz + i)] + b[getIdx(ix , iz)]) /2;
        fd_z2 += p_cur[getIdx(ix, iz-i)] * coeff[i] * (b[getIdx(ix , iz)] + b[getIdx(ix , iz - i)]) /2;
    }
    
    auto fd_sum = ((fd_x1 + fd_x2) / m_dx2) + ((fd_z1 + fd_z2) / m_dz2);

    return fd_sum;
}


/*************************************************
 * Update phi_x and phi_z
 * (1 + (zeta_x(i) + zeta_x(i+1)) / 4 * dt) * phi(t+1) = 
 * phi(t) - (zeta_x(i) + zeta_x(i+1)) / 4 * dt * phi(t) +
 * ((zeta_z(i+1) + zeta_z(i)) / 2 - (zeta_x(i+1) + zeta_x(i)) / 2) * v^2 * Du * dt
 * 
 * Du = ((p_new[i+1][j] + p_new[i+1][j+1] -
 *        p_new[i][j] - p_new[i][j+1]) / dx + 
 *       (p[i+1][j] + p[i+1][j+1] -
 *        p[i][j] - p[i][j+1]) / dx) / 4
 *************************************************/
void ForwardModeling::updatePhi(const int &ix, const int &iz)
{
    double Du_x = ((p_new[getIdx(ix + 1, iz)] + p_new[getIdx(ix + 1, iz + 1)] -
                   p_new[getIdx(ix, iz)] - p_new[getIdx(ix, iz + 1)]) /
                      m_dx +
                  (p_cur[getIdx(ix + 1, iz)] + p_cur[getIdx(ix + 1, iz + 1)] -
                   p_cur[getIdx(ix, iz)] - p_cur[getIdx(ix, iz + 1)]) /
                      m_dx) /
                 4.;
    double Du_z = ((p_new[getIdx(ix, iz + 1)] + p_new[getIdx(ix + 1, iz + 1)] -
                   p_new[getIdx(ix, iz)] - p_new[getIdx(ix + 1, iz)]) /
                      m_dz +
                  (p_cur[getIdx(ix, iz + 1)] + p_cur[getIdx(ix + 1, iz + 1)] -
                   p_cur[getIdx(ix, iz)] - p_cur[getIdx(ix + 1, iz)]) /
                      m_dz) /
                 4.;
    double term1_x = -(zeta_x[getIdx(ix, iz)] + zeta_x[getIdx(ix + 1, iz)]) / 4. * m_dt * phi_cur_x[getIdx(ix, iz)];
    double term1_z = -(zeta_z[getIdx(ix, iz)] + zeta_z[getIdx(ix, iz + 1)]) / 4. * m_dt * phi_cur_z[getIdx(ix, iz)];
    double term2_x = ((zeta_z[getIdx(ix, iz)] + zeta_z[getIdx(ix, iz + 1)]) / 2. -
                     (zeta_x[getIdx(ix, iz)] + zeta_x[getIdx(ix + 1, iz)]) / 2.) *
                    Du_x * m_dt * v[getIdx(ix, iz)] * v[getIdx(ix, iz)];
    double term2_z = ((zeta_x[getIdx(ix, iz)] + zeta_x[getIdx(ix + 1, iz)]) / 2. -
                     (zeta_z[getIdx(ix, iz)] + zeta_z[getIdx(ix, iz + 1)]) / 2.) *
                    Du_z * m_dt * v[getIdx(ix, iz)] * v[getIdx(ix, iz)];
    phi_new_x[getIdx(ix, iz)] = (phi_cur_x[getIdx(ix, iz)] + term1_x + term2_x) /
                                (1 + (zeta_x[getIdx(ix, iz)] + zeta_x[getIdx(ix + 1, iz)]) / 4. * m_dt);
    phi_new_z[getIdx(ix, iz)] = (phi_cur_z[getIdx(ix, iz)] + term1_z + term2_z) /
                                (1 + (zeta_z[getIdx(ix, iz)] + zeta_z[getIdx(ix, iz + 1)]) / 4. * m_dt);
}

/*************************************************
 * Propagate the wavefield for one step
**************************************************/
std::vector<double> ForwardModeling::step(double t)
{
    // compute the force term
    if(m_back_fd){
        setbackForce(t);
    }else{
        setForce(t);
    }
    //setForce(t);
    int idxbc = 0;
    //int idxbc2 = 0;
    //int idxbc3 = 0;
    for (int ix{half_fd_order}; ix < m_global_nx - half_fd_order; ix++)
        for (int iz{half_fd_order}; iz < m_global_nz - half_fd_order; iz++)
        {
            double idx = getIdx(ix, iz);
            double fd_sum = laplacian(ix, iz);
            double term1 = 0.;
            if(m_back_fd){
                term1 = v[idx] * v[idx] * fd_sum;
            }else{
                //term1 = v[idx] * v[idx] * (force[idx] + fd_sum);            
                double fd_sum = laplacian2(ix, iz);
                term1 = a[idx] * b[idx] * force[idx] + a[idx] * fd_sum;
            }
            
            if (isBoundary(ix, iz))
            {
                // phi[i][j] = phi_{i+1/2, j+1/2}
                // (1 + (zeta_x + zeta_z) / 2 * dt) p(t+1) = 2 * p(t) - p(t-1) +
                //      dt^2 * (v^2 * (laplacian + f) +
                //              (zeta_x + zeta_z) * p(t-1) / 2 / dt -
                //              zeta_x * zeta_z * p(t) +
                //              divergence)
                // divergence = ((phi_x[i  ][j] + phi_x[i][j-1]) / 2 - (phi_x[i-1][j  ] + phi_x[i-1][j-1]) / 2) / dx +
                //              ((phi_z[i-1][j] + phi_z[i][j  ]) / 2 - (phi_z[i-1][j-1] + phi_z[i  ][j-1]) / 2) / dz
                double term2 = (zeta_x[idx] + zeta_z[idx]) * p_old[idx] / 2. / m_dt;
                double term3 = -zeta_x[idx] * zeta_z[idx] * p_cur[idx];
                double term4 = ((phi_cur_x[getIdx(ix, iz)] + phi_cur_x[getIdx(ix, iz - 1)]) / 2 -
                               (phi_cur_x[getIdx(ix - 1, iz)] + phi_cur_x[getIdx(ix - 1, iz - 1)]) / 2) /
                                  m_dx +
                              ((phi_cur_z[getIdx(ix, iz)] + phi_cur_z[getIdx(ix - 1, iz)]) / 2 -
                               (phi_cur_z[getIdx(ix, iz - 1)] + phi_cur_z[getIdx(ix - 1, iz - 1)]) / 2) /
                                  m_dz;
                p_new[idx] = (2. * p_cur[idx] - p_old[idx] + m_dt2 * (term1 + term2 + term3 + term4)) /
                             (1. + (zeta_x[idx] + zeta_z[idx]) / 2. * m_dt);
            }else if(isTopBoundary(ix,iz)){
                    int tmp_idx = getIdx(ix+1, iz);
                    if(m_back_fd){
                        p_new[idx] = force[idx];
                    }else{
                        p_new[idx] = p_cur[tmp_idx];
                    }
                    //if( t == 0.799){
                    //    printf("force[%d, %d] = %f\n", ix, iz, force[idx]);
                    //}
                    //bcp_new[idxbc]  = p_new[idx] + g_random_nosie(m_mean, m_stddev);
                    bcp_new[idxbc]  = p_new[idx];
                    idxbc += 1;
            }
            else
                // p(t+1) = 2 * p(t) - p(t-1) + dt^2 * v^2 * (laplacian + f)
                p_new[idx] = 2. * p_cur[idx] - p_old[idx] + term1 * m_dt2;

        }

    // update Phi
    // left boundary
    for (int ix{half_fd_order}; ix < start_x; ix++)
        for (int iz{half_fd_order}; iz < m_global_nz - half_fd_order; iz++)
            updatePhi(ix, iz);
    // right boundary
    for (int ix{m_global_nx - start_x}; ix < m_global_nx - half_fd_order; ix++)
        for (int iz{half_fd_order}; iz < m_global_nz - half_fd_order; iz++)
            updatePhi(ix, iz);
    // top boundary
    for (int ix{half_fd_order}; ix < m_global_nx - half_fd_order; ix++)
        for (int iz{half_fd_order}; iz < start_z; iz++)
            updatePhi(ix, iz);
    // left boundary
    for (int ix{half_fd_order}; ix < m_global_nx - half_fd_order; ix++)
        for (int iz{m_global_nz - start_z}; iz < m_global_nz - half_fd_order; iz++)
            updatePhi(ix, iz);

    // update wavefield
    p_old = p_cur;
    p_cur = p_new;
    // update phi
    phi_cur_x = phi_new_x;
    phi_cur_z = phi_new_z;

    // write the wavefield
    if (save){
        //std::ofstream wavefieldf{"Fullwavefield_ab_2.dat", std::ios::binary | std::ios::app};
        //wavefieldf.write((char *)&p_new[0], p_new.size() * sizeof(double));
        //write the boundary wavedata
        if(m_back_fd){
            //std::ofstream wavefieldfo{"Fullwavefield_back_old.dat", std::ios::binary | std::ios::app};
            //wavefieldfo.write((char *)&p_old[0],p_old.size() * sizeof(double));

            //std::ofstream wavefieldfc{"Fullwavefield_back_cur.dat", std::ios::binary | std::ios::app};
            //wavefieldfc.write((char *)&p_cur[0],p_cur.size() * sizeof(double));

            std::ofstream wavefieldf{"Fullwavefield_back_mute.dat", std::ios::binary | std::ios::app};
            wavefieldf.write((char *)&p_new[0],p_new.size() * sizeof(double));

            std::ofstream BCwavefieldf{"Boundaryfield_back_mute.dat", std::ios::binary | std::ios::app};
            BCwavefieldf.write((char *)&bcp_new[0],bcp_new.size() * sizeof(double));
        }else{
            std::ofstream wavefieldf{"Fullwavefield.dat", std::ios::binary | std::ios::app};
            wavefieldf.write((char *)&p_new[0],p_new.size() * sizeof(double));

            std::ofstream BCwavefieldf{"Boundaryfield.dat", std::ios::binary | std::ios::app};
            BCwavefieldf.write((char *)&bcp_new[0],bcp_new.size() * sizeof(double));
        }
    }
    return bcp_new;
}

void ForwardModeling::runstep(){
    for (int it{0}; it < m_nt; it++)
    {
        double t = it * m_dt;
        bcstep = step(t);
        bcstep.clear();
    }
}
