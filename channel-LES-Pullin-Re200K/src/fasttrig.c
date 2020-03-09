#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>


#define K_MSB_BITS 15
#define K_LSB_BITS 15
#define K_MSB_TABLE_SIZE (1 << K_MSB_BITS)
#define K_LSB_TABLE_SIZE (1 << K_LSB_BITS)
#define K_BITS  (K_MSB_BITS + K_LSB_BITS)


struct tsincos
{
    double m_sin, m_cos;
};

static struct tsincos msb_table[K_MSB_TABLE_SIZE];
static struct tsincos lsb_table[K_LSB_TABLE_SIZE];
static int initialized = 0;

void fasttrig_init()
{
    int i;
    for (i = 0; i < K_MSB_TABLE_SIZE; i++) {
        double theta = (2 * M_PI * i / K_MSB_TABLE_SIZE);
        msb_table[i].m_sin = sin(theta);
        msb_table[i].m_cos = cos(theta);
    }
    
    for (i = 0; i < K_LSB_TABLE_SIZE; i++) {
        double theta = (2 * M_PI * i / (K_MSB_TABLE_SIZE * K_LSB_TABLE_SIZE));
        lsb_table[i].m_sin = sin(theta);
        lsb_table[i].m_cos = cos(theta);
    }
    printf("fasttrig_init() LUT initialized\n");
    initialized = 1;
}

/** compute sincos, accurate to one in 1<<(K_MSB_BITS + K_LSB_BITS - 3) **/
void fasttrig_sincos(double theta, double *s, double *c)
{
    if (!initialized)
        fasttrig_init();
    
    uint32_t idx = (K_MSB_TABLE_SIZE * K_LSB_TABLE_SIZE / (2*M_PI)) * theta;
    
    // rewrite theta = M + L, where L is very small.
    int lsb_idx = idx & (K_LSB_TABLE_SIZE - 1);
    double sinL, cosL;
    
    // compute sinL/cosL using lookup table
    sinL = lsb_table[lsb_idx].m_sin;
    cosL = lsb_table[lsb_idx].m_cos;
    
    int msb_idx = (idx >> K_LSB_BITS) & (K_MSB_TABLE_SIZE - 1);
    double sinM = msb_table[msb_idx].m_sin;
    double cosM = msb_table[msb_idx].m_cos;
    
    // angle sum formulas
    // we lose a few bits of precision here... about 3
    *s = sinL*cosM + sinM*cosL;
    *c = cosL*cosM - sinL*sinM;
}


/** compute sin, accurate to one in 1<<(K_MSB_BITS + K_LSB_BITS - 3) **/
double fasttrig_sin(double theta)
{
    if (!initialized)
        fasttrig_init();
    
    uint32_t idx = (K_MSB_TABLE_SIZE * K_LSB_TABLE_SIZE / (2*M_PI)) * theta;
    
    // rewrite theta = M + L, where L is very small.
    int lsb_idx = idx & (K_LSB_TABLE_SIZE - 1);
    double sinL, cosL;
    
    // compute sinL/cosL using lookup table
    sinL = lsb_table[lsb_idx].m_sin;
    cosL = lsb_table[lsb_idx].m_cos;
    
    int msb_idx = (idx >> K_LSB_BITS) & (K_MSB_TABLE_SIZE - 1);
    double sinM = msb_table[msb_idx].m_sin;
    double cosM = msb_table[msb_idx].m_cos;
    
    // angle sum formulas
    // we lose a few bits of precision here... about 3
    return sinL*cosM + sinM*cosL;
}


double fasttrig_cos(double theta)
{
    if (!initialized)
        fasttrig_init();
    
    uint32_t idx = (K_MSB_TABLE_SIZE * K_LSB_TABLE_SIZE / (2*M_PI)) * theta;
    
    // rewrite theta = M + L, where L is very small.
    int lsb_idx = idx & (K_LSB_TABLE_SIZE - 1);
    double sinL, cosL;
    
    // compute sinL/cosL using lookup table
    sinL = lsb_table[lsb_idx].m_sin;
    cosL = lsb_table[lsb_idx].m_cos;
    
    int msb_idx = (idx >> K_LSB_BITS) & (K_MSB_TABLE_SIZE - 1);
    double sinM = msb_table[msb_idx].m_sin;
    double cosM = msb_table[msb_idx].m_cos;
    
    // angle sum formulas
    // we lose a few bits of precision here... about 3
    return cosL*cosM - sinL*sinM;
}




void fasttrig_sincos_test()
{
    fasttrig_init();
    
    double eps = 1.0/(1<<17);
    
    int limit = 30;
    // int limit = INT32_MAX;
    
    int iters;
    for (iters = 0; iters < limit; iters++) {
        double theta = rand()/ 1000.0;
        double s,  c;
        double s2, c2;
        s = sin(theta);
        c = cos(theta);
        fasttrig_sincos(theta, &s2, &c2);
        
        if (1) {
            double s_err = fabs(s-s2), c_err = fabs(c-c2);
            if (s_err > eps)
                printf("sin %e : %e %e %e %10.5f\n", theta, s, s2, s_err/s, s_err/eps);
            if (c_err > eps)
                printf("cos %e : %e %e %e %10.5f\n", theta, c, c2, c_err/c, c_err/eps);
        }
    }
    exit(0);
    
}


/*
 int main(int argc, char *argv[])
 {
 fasttrig_sincos_test();
 }
 */
