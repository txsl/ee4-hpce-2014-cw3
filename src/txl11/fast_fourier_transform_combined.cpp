#include "fourier_transform.hpp"

#include <cmath>
#include <cassert>
#include <complex>

#include "tbb/parallel_for.h"
#include "tbb/task_group.h"

namespace hpce
{

namespace txl11
{

class fast_fourier_transform_combined
    : public fourier_transform
{

private:
    mutable unsigned chunk_size = 0;
    mutable unsigned recursion_size = 0;

protected:
    /* Standard radix-2 FFT only supports binary power lengths */
    virtual size_t calc_padded_size(size_t n) const
    {
        assert(n!=0);
        
        size_t ret=1;
        while(ret<n){
            ret<<=1;
        }
        
        return ret;
    }
    
    virtual size_t get_recursion_size(void) const
    {
        // printf("In get_recursion_size. val is %i\n", recursion_size);
        if (recursion_size != 0){
            return recursion_size;
        }

        char *v=getenv("HPCE_FFT_RECURSION_K");
        if(v==NULL){
            recursion_size = 16;
            // printf("HPCE_FFT_RECURSION_K not set. Using a size of %i instead.\n", recursion_size);
            
        }else{
            recursion_size = atoi(v);
            // printf("Using a chunk size of %i (set in the environment variable 'HPCE_FFT_RECURSION_K'.\n)", recursion_size);
        }
        return recursion_size;
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
        size_t n,   const std::complex<double> &wn,
        const std::complex<double> *pIn, size_t sIn,
        std::complex<double> *pOut, size_t sOut
    ) const 
    {
        assert(n>0);
     
        unsigned K_chunk_size = get_chunk_size();
        unsigned K_recursion_size = get_recursion_size();

        if (n == 1){
            pOut[0] = pIn[0];
        }else if (n == 2){
            pOut[0] = pIn[0]+pIn[sIn];
            pOut[sOut] = pIn[0]-pIn[sIn];
        }else{
            size_t m = n/2;

            if(n <= K_recursion_size){
                forwards_impl(m,wn*wn,pIn,2*sIn,pOut,sOut);
                forwards_impl(m,wn*wn,pIn+sIn,2*sIn,pOut+sOut*m,sOut);                
            }
            else{
                size_t m = n/2;

                tbb::task_group group;

                group.run( [&]() {forwards_impl(m,wn*wn,pIn,2*sIn,pOut,sOut);} );
                group.run( [&]() {forwards_impl(m,wn*wn,pIn+sIn,2*sIn,pOut+sOut*m,sOut);} );

                group.wait();
            }
             
            if(m <= K_chunk_size){

                std::complex<double> w=std::complex<double>(1.0, 0.0);

                for (size_t j=0;j<m;j++){
                    std::complex<double> t1 = w*pOut[m+j];
                    std::complex<double> t2 = pOut[j]-t1;
                    pOut[j] = pOut[j]+t1;                 /*  pOut[j] = pOut[j] + w^i pOut[m+j] */
                    pOut[j+m] = t2;                          /*  pOut[j] = pOut[j] - w^i pOut[m+j] */
                    w = w*wn;
                }
            }
            else{
                typedef tbb::blocked_range<unsigned> my_range_t;

                my_range_t range(0, m, K_chunk_size);

                auto f = [=](const my_range_t &chunk){
                    std::complex<double> w_local = std::complex<double>(1.0, 0.0);
                    w_local = pow(wn, chunk.begin());

                    for (size_t j=chunk.begin(); j!=chunk.end(); j++){

                        std::complex<double> t1 = w_local*pOut[m+j];
                        std::complex<double> t2 = pOut[j]-t1;
                        pOut[j] = pOut[j]+t1;                 /*  pOut[j] = pOut[j] + w^i pOut[m+j] */
                        pOut[j+m] = t2;                          /*  pOut[j] = pOut[j] - w^i pOut[m+j] */
                        w_local = w_local*wn;
                    }
                };

                tbb::parallel_for (range, f, tbb::simple_partitioner());
               
            }
        }
    }
    
    virtual void backwards_impl(
        size_t n,   const std::complex<double> &wn,
        const std::complex<double> *pIn, size_t sIn,
        std::complex<double> *pOut, size_t sOut
    ) const 
    {
        complex_t reverse_wn=1.0/wn;
        forwards_impl(n, reverse_wn, pIn, sIn, pOut, sOut);
        
        double scale=1.0/n;
        for(size_t i=0;i<n;i++){
            pOut[i]=pOut[i]*scale;
        }
    }
    
public:
    virtual std::string name() const
    { return "hpce.txl11.fast_fourier_transform_combined"; }
    
    virtual bool is_quadratic() const
    { return false; }
};

std::shared_ptr<fourier_transform> Create_fast_fourier_transform_combined()
{
    return std::make_shared<fast_fourier_transform_combined>();
}

}; // namespace txl11
}; // namespace hpce
