/***************************************************************************
 *                                                                         *
 * Arm-Rotor-Platform Example                                              *
 *                                                                         *
 * Authors:                                                                *
 *   Juergen Pannek <juergen.pannek@uni-bayreuth.de>                       *
 *                                                                         *
 ***************************************************************************/

/**
 @file
 @author Juergen Pannek <juergen.pannek@uni-bayreuth.de>
 */

#ifndef ARP_H
#define ARP_H

#include <yane/model.h>

namespace yane
{
	namespace ExampleModels2
	{
		class ARP : public yane::Model::Model
		{
			public:
				ARP ( );
				~ARP ( );

				virtual void odeFunction ( double t, double * x, double * u,
				        double * dx );
				virtual double continuousCostfunction ( double t, double * x,
				        double * u );
				virtual void restrictionFunction ( double t, double *x,
				        double * u, double * fx );
				virtual double discreteCostfunction ( int length, int horizont,
				        double * t, double * x, double * u );
				virtual void getObjectiveWeight (
				        double & weight_continuousCostfunction,
				        double & weight_discreteCostfunction );
				virtual void getDefaultState ( double * x );
				virtual void getDefaultControl ( double * u );
				virtual void eventBeforeMPC ( int horizont, double * timesteps,
				        double * x, double * sdatavalues );

				virtual const yane::Utils::Uuid & id ( );
				static const yane::Utils::Uuid ID;

			private:
				void setModelParameter ( void );
		};
	}
}

#endif
