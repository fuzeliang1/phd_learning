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

#include "../Models/model_7_31.h"

#include <yane.h>
#include <sqpfortran.h>

#define HORIZON 50
#define H_NEW 0.15
#define INF 1.0E19
#define HORIZONSHIFT 1

using namespace yane::MPC;
using namespace yane::SqpFortran;
using namespace yane::Utils;
using namespace yane::ExampleModels2;
using namespace std;

int main ( )
{
	double * x, *u, *t, *nextControls, *nextTimesteps;

	yane::Model::Model * model = new Crane ( );
	yane::Model::Simulator * simulator = new yane::Model::Simulator ( model );

	OdeManager * odemanager = new CacheOdeManager ( );
	yane::MinProg::NLP * sqp = new SqpFortran ( );

	SqpFortran * sqpSpecify = ( ( SqpFortran * ) sqp );
	// 	sqpSpecify->setAccuracy ( 1E-4 );
	// 	sqpSpecify->setMaxFun ( 20 );
	// 	sqpSpecify->setMaxIterations ( 1000 );
	// 	sqpSpecify->setLineSearchTol ( 0.1 );

	MPC * mpc = new SuboptimalityMPC ( H_NEW, INF );
	SuboptimalityMPC * mpcSpecify = ( SuboptimalityMPC * ) mpc;
	mpc->reset ( odemanager, sqp, model, HORIZON, 0 );

	mpc->allocateMemory ( t, u );
	mpc->generateTimeGrid ( t, 0.0, H_NEW );

	mpc->initCalc ( t, u );
	mpc->resizeHorizon ( 17, H_NEW );
	mpc->setConfigShiftHorizon ( STANDARD );

	OptimalDiff::calc ( );
	Exception::enableDebugMessage ( true );

	x = new double [ model->dimensionState ( ) ];

	nextControls = new double [ model->dimensionControl ( ) * HORIZONSHIFT ];
	nextTimesteps = new double [ HORIZONSHIFT + 1 ];

	int index, mpciterations;
	double closedloopcost, alphamax, alpha, gamma;
	double * valuefunction = new double [ HORIZON + 1 ];
	double * trajectory = new double [ ( HORIZON + 1 )
	        * model->dimensionState ( ) ];

	double alphabound [ ] =
	{ 0.1, 0.5, 0.9, 0.95 };
	for ( int alphaiteration = 0; alphaiteration < 4; alphaiteration++ )
	{
		mpciterations = 0;
		closedloopcost = 0.0;
		model->getDefaultState ( x );
		mpc->generateTimeGrid ( t, 0.0, H_NEW );
		for ( int i = 0; i < HORIZON; i++ )
		{
			model->getDefaultControl ( & u [ i * model->dimensionControl ( ) ] );
		}

		cout
		        << "  Time          u       phi      omega          x          v       alpha      ms"
		        << endl;
		cout
		        << "================================================================================"
		        << endl;
		while ( ( double ) mpciterations * H_NEW <= 15.0 )
		{
			RTClock timer;
			try
			{
				mpc->calc ( x, 0, 0.0, alphabound [ alphaiteration ] );
			}
			catch ( yane::MinProg::SolverWarning e )
			{
				// 				cout << e.what() << endl;
			}

			odemanager->getObjFunValues ( t, x, u, valuefunction, trajectory );
			closedloopcost += valuefunction [ 0 ];

			alphamax = 1E-20;
			try
			{
				// 				mpcSpecify->aposterioriEstimate ( x, alpha );
				mpcSpecify->setEpsilon ( 1E-5 );
				mpcSpecify->aposterioriPracticalEstimate ( x, alpha );
			}
			catch ( yane::MPC::SuboptimalityException e )
			{
				// 				cout << e.what() << endl;
				index = 1;
			}
			catch ( yane::MPC::PracticalException e )
			{
				// 				cout << e.what() << endl;
				index = 2;
			}
			catch ( yane::MinProg::SolverWarning e )
			{
				// 				cout << e.what() << endl;
			}

			double duration = timer.elapsedSeconds ( ) * 1000.0;

			try
			{
				mpc->shiftHorizon ( nextControls, nextTimesteps, HORIZONSHIFT,
				        x, H_NEW );
			}
			catch ( yane::MinProg::SolverWarning e )
			{
				// 				cout << e.what() << endl;
			}

			for ( int j = 0; j < HORIZONSHIFT; j++ )
			{
				cout << fixed << setprecision ( 2 ) << setw ( 6 )
				        << nextTimesteps [ j ];
				cout << fixed << setprecision ( 5 ) << setw ( 11 )
				        << nextControls [ j ];
				for ( int k = 0; k < 4; k++ )
				{
					cout << setprecision ( 5 ) << setw ( 10 ) << x [ k ] << " ";
				}
				cout << fixed << setprecision ( 5 ) << setw ( 11 ) << alpha;
				if ( j == 0 )
				{
					cout << setw ( 8 ) << setprecision ( 2 ) << duration;
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

				simulator->initialize ( nextTimesteps [ j ], x,
				        & nextControls [ j * model->dimensionControl ( ) ] );
				simulator->run ( nextTimesteps [ j + 1 ] );
				mpciterations++ ;
			}
		}

		cout << endl << setprecision ( 5 ) << setw ( 10 ) << closedloopcost
		        << endl << endl << endl;
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
	delete [ ] valuefunction;
	delete [ ] trajectory;

	return 0;
}
