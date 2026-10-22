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

#ifndef PENDULUMSHOOTINGDATA_H
#define PENDULUMSHOOTINGDATA_H

#include <yane/mpc/modelshootingdata.h>

namespace yane
{
	namespace Model
	{
		class Model;
	}

	namespace ExampleModels2
	{
		class InvertedPendulumShootingData : public yane::MPC::ModelShootingData
		{
			public:
				InvertedPendulumShootingData ( yane::Model::Model * pendulum );
				virtual ~InvertedPendulumShootingData ( );

				virtual void applyShootingData ( int horizon, double * t,
				        double * x, double * sdatavalues );
				virtual int shootingDataLength ( int horizon );
				virtual int maxShootingDataLength ( int maxhorizon );
				virtual void shootingDataInfo ( int horizon,
				        yane::MPC::T_SHOOTINGDATA * sdata );

				virtual void setShooting ( int horizon,
				        int numberOfShootingNodes, int * shootingHorizonIndex,
				        int * shootingVarIndex, double * shootingValue );

			protected:
				yane::Model::Model * _pendulum;
				int _numberOfShootingNodes;
				int * _shootingHorizonIndex;
				int * _shootingVarIndex;
				double * _shootingValue;
		};
	}
}

#endif
