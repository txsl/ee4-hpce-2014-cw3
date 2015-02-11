#include "fourier_transform.hpp"

namespace hpce{

// Declare factory functions which are implemented elsewhere.
std::shared_ptr<fourier_transform> Create_fast_fourier_transform();
std::shared_ptr<fourier_transform> Create_direct_fourier_transform();
	
// TODO : Declare your factories here

namespace txl11{
	std::shared_ptr<fourier_transform> Create_direct_fourier_transform_parfor();
	std::shared_ptr<fourier_transform> Create_fast_fourier_transform_taskgroup();
}



void fourier_transform::RegisterDefaultFactories()
{
	static const unsigned MYSTERIOUS_LINE=0; // Don't remove me!
	
	RegisterTransformFactory("hpce.fast_fourier_transform", Create_fast_fourier_transform);
	RegisterTransformFactory("hpce.direct_fourier_transform", Create_direct_fourier_transform);
	
	// TODO : Add your factories here
	// e.g. RegisterTransformFactory("hpce.YOUR_LOGIN.direct_fourier_transform_parfor", hpce::YOUR_LOGIN::Create_direct_fourier_transform_parfor);
	RegisterTransformFactory("hpce.txl11.direct_fourier_transform_parfor", hpce::txl11::Create_direct_fourier_transform_parfor);
	RegisterTransformFactory("hpce.txl11.fast_fourier_transform_taskgroup", hpce::txl11::Create_fast_fourier_transform_taskgroup);
}
	
}; // namespace hpce
