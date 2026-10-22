/***************************************************************************
 *                                                                         *
 * Inverted Pendulum on a Cart Example                                     *
 *                                                                         *
 * Authors:                                                                *
 *   Juergen Pannek <juergen.pannek@uni-bayreuth.de>                       *
 *                                                                         *
 ***************************************************************************/

/**
 @file
 @author Juergen Pannek <juergen.pannek@uni-bayreuth.de>
 */

#ifndef INVERTEDPENDULUM2_H
#define INVERTEDPENDULUM2_H

#include <yane/model.h>

namespace yane
{
	namespace ExampleModels2
	{

		class InvertedPendulum2 : public yane::Model::Model
		{
			public:
				InvertedPendulum2 ( );
				~InvertedPendulum2 ( );

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

				double getTerminalWeight ( void );
				void setTerminalWeight ( double terminal_weight,
				        double * terminal_time );

			private:
				double _terminalWeight;
				double * _terminalTime;

		};
	}
}

#endif
