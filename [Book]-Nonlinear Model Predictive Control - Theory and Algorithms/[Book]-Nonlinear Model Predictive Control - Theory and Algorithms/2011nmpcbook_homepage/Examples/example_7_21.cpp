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
#include <fstream>

#include "../Models/model_7_21.h"

#include <yane.h>
#include <sqpfortran.h>

#define HORIZON 20
#define H_NEW 0.2
#define INF 1.0E19
#define HORIZONSHIFT 1

using namespace yane::MPC;
using namespace yane::SqpFortran;
using namespace yane::Utils;
using namespace yane::ExampleModels2;
using namespace yane::Model;
using namespace std;

int main ( )
{
	double * x, *u, *t, *nextControls, *nextTimesteps;

	yane::Model::Model * model = new ARP ( );
	yane::Model::Simulator * simulator = new yane::Model::Simulator ( model );

	OdeManager * odemanager = new CacheOdeManager ( );
	yane::MinProg::NLP * sqp = new SqpFortran ( );

	SqpFortran * sqpSpecify = ( ( SqpFortran * ) sqp );
	sqpSpecify->setAccuracy ( 1E-6 );
	sqpSpecify->setMaxFun ( 20 );
	sqpSpecify->setMaxIterations ( 1000 );
	sqpSpecify->setLineSearchTol ( 0.1 );

	MPC * mpc = new AdaptiveMPC ( H_NEW, INF );
	AdaptiveMPC * mpcSpecify = ( ( AdaptiveMPC * ) mpc );
	mpc->reset ( odemanager, sqp, model, HORIZON, 0 );

	double alphabound = 0.1;
	int N0 = 2;
	int sigma = 5;
	double eps = 1E-3;
	mpcSpecify->setConfigAdaptiveMPC ( alphabound, APOSTERIORIPRACTICAL,
	        SIMPLE, N0, eps, sigma, true, true );
	int N = 7;

	mpc->allocateMemory ( t, u );
	mpc->generateTimeGrid ( t, 0.0, H_NEW );
	for ( int i = 0; i < HORIZON; i++ )
	{
		model->getDefaultControl ( & u [ i * model->dimensionControl ( ) ] );
	}

	mpc->initCalc ( t, u );
	mpc->resizeHorizon ( N, H_NEW );
	mpc->setConfigShiftHorizon ( STANDARD );

	nextControls = new double [ model->dimensionControl ( ) * HORIZONSHIFT ];
	nextTimesteps = new double [ HORIZONSHIFT + 1 ];
	x = new double [ model->dimensionState ( ) ];
	model->getDefaultState ( x );

	OptimalDiff::calc ( );
	Exception::enableDebugMessage ( true );

	int mpciterations = 0;
	double computingtime = 0.0, totalcomputingtime = 0.0;

	cout
	        << "    Zeit          u       x(9)      x(10)      alpha        ms   N"
	        << endl;
	cout
	        << "=================================================================="
	        << endl;
	while ( ( double ) mpciterations * H_NEW <= 15.0 )
	{
		RTClock timer;
		try
		{
			mpcSpecify->calc ( x );
		}
		catch ( yane::MinProg::SolverWarning e )
		{
			cout << e.what ( ) << endl;
		}
		catch ( yane::Utils::Exception e )
		{
			cout << e.what ( ) << endl;
		}

		computingtime = timer.elapsedSeconds ( ) * 1000.0;
		totalcomputingtime += computingtime;

		try
		{
			mpcSpecify->shiftHorizon ( nextControls, nextTimesteps,
			        HORIZONSHIFT, x, H_NEW );
		}
		catch ( yane::MinProg::SolverWarning e )
		{
			cout << e.what ( ) << endl;
		}

		for ( int j = 0; j < HORIZONSHIFT; j++ )
		{
			cout << fixed << setprecision ( 2 ) << setw ( 8 )
			        << nextTimesteps [ 0 ];
			cout << fixed << setprecision ( 5 ) << setw ( 11 )
			        << nextControls [ 0 ];
			cout << fixed << setprecision ( 5 ) << setw ( 11 ) << x [ 8 ];
			cout << fixed << setprecision ( 5 ) << setw ( 11 ) << x [ 9 ];
			cout << fixed << setprecision ( 5 ) << setw ( 11 )
			        << mpcSpecify->alpha ( );
			if ( j == 0 )
			{
				cout << setw ( 10 ) << setprecision ( 2 ) << computingtime;
			}
			else
			{
				cout << "        --";
			}
			cout << fixed << setprecision ( 0 ) << setw ( 4 )
			        << mpc->horizon ( );
			cout << endl;

			simulator->initialize ( nextTimesteps [ j ], x, & nextControls [ j
			        * model->dimensionControl ( ) ] );
			simulator->run ( nextTimesteps [ j + 1 ] );
			mpciterations++ ;
		}
	}
	cout << setprecision ( 2 ) << setw ( 10 ) << endl << totalcomputingtime
	        << endl;

	delete odemanager;
	delete sqp;
	delete model;
	delete mpc;
	delete [ ] x;
	delete [ ] u;
	delete [ ] t;
	delete [ ] nextControls;
	delete [ ] nextTimesteps;

	return 0;

}

