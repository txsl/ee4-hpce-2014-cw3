#include "fourier_transform.hpp"

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <string>

#include "tbb/parallel_for.h"

namespace hpce
{

namespace txl11 {

class direct_fourier_transform_chunked
    : public fourier_transform
{
private:
    mutable unsigned chunk_size = 0;

protected:
    /*! We can do any size transform */
    virtual size_t calc_padded_size(size_t n) const
    {
        return n;
    }

    virtual size_t get_chunk_size(void) const
    {
        // printf("In get_recursion_size. val is %i\n",chunk_size);
        if (chunk_size != 0){
            return chunk_size;
        }

        char *v=getenv("HPCE_FFT_LOOP_K");
        if(v==NULL){
           chunk_size = 16;
            // printf("HPCE_FFT_LOOP_K not set. Using a size of %i instead.\n", chunk_size);
            
        }else{
           chunk_size = atoi(v);
            // printf("Using a chunk size of %i (set in the environment variable 'HPCE_FFT_LOOP_K'.\n)", chunk_size);
        }
        return chunk_size;
    }

    virtual void forwards_impl(
        size_t n,   const std::complex<double> &/*wn*/,
        const std::complex<double> *pIn, size_t sIn,
        std::complex<double> *pOut, size_t sOut
    ) const 
    {

        assert(n>0);
        
        const double PI=3.1415926535897932384626433832795;
        
        // = -i*2*pi / n
        complex_t neg_im_2pi_n=-complex_t(0.0, 1.0)*2.0*PI / (double)n;

        unsigned K = get_chunk_size();

        // As used in an example in the coursework spec
        typedef tbb::blocked_range<unsigned> my_range_t;

        my_range_t range(0, n, K);

        auto f = [=](const my_range_t &chunk){
            for(unsigned kk=chunk.begin(); kk!=chunk.end(); kk++ ){
                complex_t acc=0;

                for(size_t ii=0;ii<n;ii++){
                    // acc += exp(-i * 2 * pi * kk / n);
                    acc+=pIn[ii*sIn] * exp( neg_im_2pi_n * (double)kk * (double)ii );
                }

                pOut[kk*sOut]=acc;
                }
        };

        tbb::parallel_for (range, f, tbb::simple_partitioner());
    }
    
    virtual void backwards_impl(
        size_t n,   const std::complex<double> &/*wn*/,
        const std::complex<double> *pIn, size_t sIn,
        std::complex<double> *pOut, size_t sOut
    ) const 
    {
        assert(n>0);
        
        // unsigned K = get_chunk_size();

        const double PI=3.1415926535897932384626433832795;
        
        // = i*2*pi / n
        complex_t im_2pi_n=complex_t(0.0, 1.0)*2.0*PI / (double)n;
        
        const double scale=1.0/n;

        unsigned K = get_chunk_size();

        typedef tbb::blocked_range<unsigned> my_range_t;

        my_range_t range(0, n, K);

        auto f = [=](const my_range_t &chunk){
            for(unsigned kk=chunk.begin(); kk!=chunk.end(); kk++ ){
                complex_t acc=0;

                for(size_t ii=0;ii<n;ii++){
                    // acc += exp(-i * 2 * pi * kk / n);
                    acc+=pIn[ii*sIn] * exp( im_2pi_n * (double)kk * (double)ii );
                }

                pOut[kk*sOut]=acc*scale;
                }
        };

        tbb::parallel_for (range, f, tbb::simple_partitioner());
    }
    
public:
    virtual std::string name() const
    { return "hpce.txl11.direct_fourier_transform_chunked"; }
    
    virtual bool is_quadratic() const
    { return true; }
};

std::shared_ptr<fourier_transform> Create_direct_fourier_transform_chunked()
{
    return std::make_shared<direct_fourier_transform_chunked>();
}

}; // namespace txl11
}; // namespace hpce
