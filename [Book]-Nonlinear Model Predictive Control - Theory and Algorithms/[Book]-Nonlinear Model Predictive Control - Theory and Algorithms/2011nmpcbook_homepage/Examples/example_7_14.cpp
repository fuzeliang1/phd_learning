/***************************************************************************
 *                                                                         *
 *                                 NMPC                                    *
 *                         www.nonlinearmpc.com                            *
 *                                                                         *
 * Author:                                                                *
 *   Juergen Pannek <juergen.pannek@googlemail.com>                        *
 *                                                                         *
 ***************************************************************************/

/**
 @file
 @author Juergen Pannek <juergen.pannek@googlemail.com>
 */

#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdio>
#include <unistd.h>
#include <cstdlib>

#include "../Models/model_7_14.h"

#include <yane.h>
#include <sqpfortran.h>

#define HORIZON 70
#define H_NEW 0.05
#define INF 1.0E19
#define HORIZONSHIFT 10

using namespace yane::MPC;
using namespace yane::SqpFortran;
using namespace yane::Utils;
using namespace yane::ExampleModels2;
using namespace std;

int main ( )
{
	double * x, *u, *t, *nextControls, *nextTimesteps;

	yane::Model::Model * model = new InvertedPendulum2 ( );
	InvertedPendulum2 * modelspecify = ( ( InvertedPendulum2 * ) model );
	yane::Model::Simulator * simulator = new yane::Model::Simulator ( model );

	OdeManager * odemanager = new CacheOdeManager ( );
	yane::MinProg::NLP * sqp = new SqpFortran ( );

	SqpFortran * sqpSpecify = ( ( SqpFortran * ) sqp );
	sqpSpecify->setAccuracy ( 1E-6 );
	sqpSpecify->setMaxFun ( 20 );
	sqpSpecify->setMaxIterations ( 1000 );
	sqpSpecify->setLineSearchTol ( 0.1 );

	MPC * mpc = new SuboptimalityMPC ( H_NEW, INF );
	SuboptimalityMPC * mpcSpecify = ( SuboptimalityMPC * ) mpc;
	mpc->reset ( odemanager, sqp, model, HORIZON, 0 );

	mpc->allocateMemory ( t, u );
	mpc->generateTimeGrid ( t, 0.0, H_NEW );

	mpc->initCalc ( t, u );
	mpc->setConfigShiftHorizon ( STANDARD );

	OptimalDiff::calc ( );
	Exception::enableDebugMessage ( true );

	x = new double [ model->dimensionState ( ) ];

	nextControls = new double [ model->dimensionControl ( ) * HORIZONSHIFT ];
	nextTimesteps = new double [ HORIZONSHIFT + 1 ];

	int index, mpciterations;
	double alpha;

	for ( int OptimizationHorizon = 70; OptimizationHorizon <= HORIZON; OptimizationHorizon++ )
	{
		mpc->resizeHorizon ( OptimizationHorizon, H_NEW );
		mpc->reset ( odemanager, sqp, model, OptimizationHorizon );
		for ( int ControlHorizon = 5; ControlHorizon <= HORIZONSHIFT; ControlHorizon++ )
		{
			for ( int TerminalWeight = 5; TerminalWeight <= 10; TerminalWeight++ )
			{
				mpciterations = 0;
				model->getDefaultState ( x );
				mpc->generateTimeGrid ( t, 0.0, H_NEW );
				for ( int i = 0; i < OptimizationHorizon; i++ )
				{
					model->getDefaultControl ( & u [ i
					        * model->dimensionControl ( ) ] );
				}
				mpc->initCalc ( t, u );

				cout
				        << "  Time          u       phi      omega          x          v       alpha      ms"
				        << endl;
				cout
				        << "================================================================================"
				        << endl;
				while ( ( double ) mpciterations * H_NEW <= 10.0 )
				{
					modelspecify->setTerminalWeight (
					        ( double ) TerminalWeight,
					        & t [ OptimizationHorizon - 1 ] );
					RTClock timer;
					try
					{
						mpc->calc ( x );
					}
					catch ( yane::MinProg::SolverWarning e )
					{
						cout << e.what ( ) << endl;
					}

					alpha = 1E-20;
					try
					{
						// 						mpcSpecify->mstepAposterioriEstimate ( x, alpha, ControlHorizon );
						mpcSpecify->setEpsilon ( 1E-5 );
						mpcSpecify->mstepAposterioriPracticalEstimate ( x,
						        alpha, ControlHorizon );
					}
					catch ( yane::MPC::SuboptimalityException e )
					{
						cout << e.what ( ) << endl;
						index = 1;
					}
					catch ( yane::MPC::PracticalException e )
					{
						cout << e.what ( ) << endl;
						index = 2;
					}
					catch ( yane::MinProg::SolverWarning e )
					{
						cout << e.what ( ) << endl;
					}

					double duration = timer.elapsedSeconds ( ) * 1000.0;

					try
					{
						mpc->shiftHorizon ( nextControls, nextTimesteps,
						        ControlHorizon, x, H_NEW );
					}
					catch ( yane::MinProg::SolverWarning e )
					{
						cout << e.what ( ) << endl;
					}

					for ( int j = 0; j < ControlHorizon; j++ )
					{
						cout << fixed << setprecision ( 2 ) << setw ( 6 )
						        << nextTimesteps [ j ];
						cout << fixed << setprecision ( 5 ) << setw ( 11 )
						        << nextControls [ j ];
						for ( int k = 0; k < 4; k++ )
						{
							cout << setprecision ( 5 ) << setw ( 10 ) << x [ k ]
							        << " ";
						}
						cout << fixed << setprecision ( 5 ) << setw ( 11 )
						        << alpha;
						if ( j == 0 )
						{
							cout << setw ( 8 ) << setprecision ( 2 )
							        << duration;
						}
						else
						{
							cout << "     --";
						}
						if ( index == 1 )
						{
							cout << " S";
						}
						else if ( index == 2 )
						{
							cout << " P";
						}
						cout << endl;
						index = 0;

						simulator->initialize (
						        nextTimesteps [ j ],
						        x,
						        & nextControls [ j * model->dimensionControl ( ) ] );
						simulator->run ( nextTimesteps [ j + 1 ] );
						mpciterations++ ;
					}
				}
			}
		}
	}

	delete mpc;
	delete sqp;
	delete odemanager;
	delete model;
	delete [ ] x;
	delete [ ] u;
	delete [ ] t;
	delete [ ] nextControls;
	delete [ ] nextTimesteps;

	return 0;
}

