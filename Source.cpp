/* Source for wave propagation */

#include <cmath>
#include <fstream>
#include <vector>
#include "Source.h"

const float pi{3.1415926};
Source::Source(int x, int z, float A, float f0, float t0)
    : m_x{x}, m_z{z}, m_A{A}, m_f0{f0}, m_t0{t0}
{
}

/*************************************************
 * Compute the ricker wavelet at a given time
 * f(t) = A * (1 - 2 * pi^2 * f0^2 * (t - t0)^2) 
 *      * exp(-pi^2 * f0^2 * (t - t0)^2)
 * @param t the given time to compute the ricker wavelet.
 * @return the computed amplitude.
 *************************************************/
float Source::rickerWavelet(const float &t)
{
    float x{static_cast<float>(pow(pi * m_f0 * (t - m_t0), 2))};
    return m_A * exp(-x) * (1 - 2. * x);
}

/*************************************************
 * Write the ricker wavelet to the specified file
 * @param filename the given filename to write the wavelet.
 * @param nt the length of wavelet.
 * @param dt the time interval.
**************************************************/
void Source::writeWavelet(const std::string &filename,
                          const int &nt, const float &dt)
{
    std::vector<float> wavelet{};
    std::ofstream waveletf{filename, std::ios::binary};

    for (int i{0}; i < nt; i++)
        wavelet.push_back(rickerWavelet(i * dt));

    waveletf.write((char *)&wavelet[0], wavelet.size() * sizeof(float));
}

int Source::getX()
{
    return m_x;
}

int Source::getZ()
{
    return m_z;
}
