#include "fourier_transform.hpp"

#include <cmath>
#include <cassert>

#include "tbb/task_group.h"

namespace hpce
{

namespace txl11 {

class fast_fourier_transform_taskgroup
    : public fourier_transform
{
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

    virtual void forwards_impl(
        size_t n,   const std::complex<double> &wn,
        const std::complex<double> *pIn, size_t sIn,
        std::complex<double> *pOut, size_t sOut
    ) const 
    {
        assert(n>0);
        
        if (n == 1){
            pOut[0] = pIn[0];
        }else if (n == 2){
            pOut[0] = pIn[0]+pIn[sIn];
            pOut[sOut] = pIn[0]-pIn[sIn];
        }else{
            size_t m = n/2;

            tbb::task_group group;

            group.run( [&]() {forwards_impl(m,wn*wn,pIn,2*sIn,pOut,sOut);} );
            group.run( [&]() {forwards_impl(m,wn*wn,pIn+sIn,2*sIn,pOut+sOut*m,sOut);} );

            group.wait();
             
            std::complex<double> w=std::complex<double>(1.0, 0.0);

            for (size_t j=0;j<m;j++){
              std::complex<double> t1 = w*pOut[m+j];
              std::complex<double> t2 = pOut[j]-t1;
              pOut[j] = pOut[j]+t1;                 /*  pOut[j] = pOut[j] + w^i pOut[m+j] */
              pOut[j+m] = t2;                          /*  pOut[j] = pOut[j] - w^i pOut[m+j] */
              w = w*wn;
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
    { return "hpce.txl11.fast_fourier_transform_taskgroup"; }
    
    virtual bool is_quadratic() const
    { return false; }
};

std::shared_ptr<fourier_transform> Create_fast_fourier_transform_taskgroup()
{
    return std::make_shared<fast_fourier_transform_taskgroup>();
}

}; // namespace txl11

}; // namespace hpce
