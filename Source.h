#ifndef SOURCE_H
#define SOURCE_H

class Source
{
private:
    int m_x, m_z; // location of the source
    float m_A;    // amplitude
    float m_f0;   // peak frequency
    float m_t0;   // time lag

public:
    Source(int x, int z, float A, float f0, float t0);
    float rickerWavelet(const float &t);
    void writeWavelet(const std::string &filename, 
                      const int &nt, const float &dt);
    int getX();
    int getZ();
};

#endif