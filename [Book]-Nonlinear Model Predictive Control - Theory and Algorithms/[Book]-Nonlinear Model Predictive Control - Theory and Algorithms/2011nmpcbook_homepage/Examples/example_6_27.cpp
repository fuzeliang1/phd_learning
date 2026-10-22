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

#include "../Models/model_6_27.h"

#include <yane.h>
#include <sqpfortran.h>

#define HORIZON 11
#define H_NEW 0.025
#define INF 1.0E19
#define HORIZONSHIFT 1
#define DIMENSIONSTATE 50
#define DIMENSIONCONTROL 50

using namespace yane::MPC;
using namespace yane::Model;
using namespace yane::ExampleModels2;
using namespace yane::SqpFortran;
using namespace yane::Utils;
using namespace std;

int main ( )
{
	OptimalDiff::calc ( );
	Exception::enableDebugMessage ( true );

	double * x, *u, *t;

	Model * model = new SemiHeat ( DIMENSIONSTATE, DIMENSIONCONTROL );
	Simulator * simulator = new Simulator ( model );

	OdeManager * odemanager = new MultiThreadCacheOdeManager ( 4 );
	yane::MinProg::NLP * sqp = new SqpFortran ( );

	SqpFortran * sqpSpecify = ( ( SqpFortran * ) sqp );
	sqpSpecify->setAccuracy ( 1.E-9 );
	sqpSpecify->setMaxFun ( 50 );
	sqpSpecify->setMaxIterations ( 10000 );
	sqpSpecify->setLineSearchTol ( 0.1 );

	MPC * mpc = new MPC ( INF );
	mpc->reset ( odemanager, sqp, model, HORIZON );
	mpc->setConfigShiftHorizon ( STANDARD );

	mpc->allocateMemory ( t, u );
	mpc->generateTimeGrid ( t, 0.0, H_NEW );
	mpc->initCalc ( t, u );
	x = new double [ model->dimensionState ( ) ];
	model->getDefaultState ( x );

	double * nextControls = new double [ model->dimensionControl ( )
	        * HORIZONSHIFT ];
	double * nextTimesteps = new double [ HORIZONSHIFT + 1 ];

	RTClock timer;
	cout << "  Time   u[DIM/2]          V  comp. time" << endl;
	cout << "========================================" << endl;
	for ( int i = 0; i < 41; i++ )
	{
		RTClock timer;
		try
		{
			mpc->calc ( x );
		}
		catch ( yane::MinProg::SolverWarning e )
		{
			cout << e.what ( ) << endl;
		}
		double duration = timer.elapsedSeconds ( ) * 1000.0;
		try
		{
			mpc->shiftHorizon ( nextControls, nextTimesteps, HORIZONSHIFT, x,
			        H_NEW );
		}
		catch ( yane::MinProg::SolverWarning e )
		{
			cout << e.what ( ) << endl;
		}

		for ( int j = 0; j < HORIZONSHIFT; j++ )
		{
			cout << fixed << setprecision ( 3 ) << setw ( 6 )
			        << nextTimesteps [ 0 ] << " " << setprecision ( 5 )
			        << setw ( 10 ) << nextControls [ DIMENSIONCONTROL / 2 ]
			        << " " << setprecision ( 5 ) << setw ( 10 )
			        << ( model->continuousCostfunction ( 0, x, nextControls ) );
			if ( j == 0 )
			{
				cout << setw ( 12 ) << setprecision ( 2 ) << duration;
			}
			else
			{
				cout << "        --";
			}
			cout << endl;

			simulator->initialize ( nextTimesteps [ j ], x, & nextControls [ j ] );
			simulator->run ( nextTimesteps [ j + 1 ] );
		}
	}

	delete mpc;
	delete sqp;
	delete odemanager;
	delete simulator;
	delete model;
	delete [ ] t;
	delete [ ] u;
	delete [ ] x;
	delete [ ] nextControls;
	delete [ ] nextTimesteps;

	return 0;
}
