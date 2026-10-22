/***************************************************************************
 *                                                                         *
 * Advection-Diffusion PDE                                                 *
 *                                                                         *
 * Authors:                                                                *
 *   Juergen Pannek <juergen.pannek@uni-bayreuth.de>                       *
 *                                                                         *
 ***************************************************************************/

/**
 @file
 @author Juergen Pannek <juergen.pannek@uni-bayreuth.de>
 */

#ifndef ADVECTIONDIFFUSION_H
#define ADVECTIONDIFFUSION_H

#include <yane/model.h>

namespace yane
{
	namespace ExampleModels2
	{
		class SemiHeat : public yane::Model::Model
		{
			public:
				SemiHeat ( int dimensionState, int dimensionControl );
				~SemiHeat ( );

				virtual void odeFunction ( double t, double * x, double * u,
				        double * dx );
				virtual double continuousCostfunction ( double t, double * x,
				        double * u );
				virtual void restrictionFunction ( double t, double *x,
				        double * u, double * fx );
				virtual double discreteCostfunction ( int length, int horizon,
				        double * t, double * x, double * u );
				virtual void getObjectiveWeight (
				        double & weight_continuousCostfunction,
				        double & weight_discreteCostfunction );
				virtual void getDefaultState ( double * x );
				virtual void getDefaultControl ( double * u );
				virtual void eventBeforeMPC ( int horizon, double * timesteps,
				        double * x, double * sdatavalues );

				double semiFunction ( double x );

			private:
				bool _distributedControl;
				bool _h1norm;
		};
	}
}

#endif
