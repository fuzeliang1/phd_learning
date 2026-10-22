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

#include "../Models/model_10.h"
#include "../Models/model_10_shootingdata.h"

#include <yane.h>
#include <sqpfortran.h>

#define HORIZONMAX 200
#define H_NEW 0.1
#define INF 1.0E19
#define HORIZONSHIFT 1
#define PRECEDESTEPS 0

using namespace yane::MPC;
using namespace yane::Model;
using namespace yane::ExampleModels2;
using namespace yane::SqpFortran;
using namespace yane::Utils;
using namespace std;

int main()
{
	OptimalDiff::calc();
	Exception::enableDebugMessage ( true );

	int optimizationHorizon = 20;

	double * x, * u, * t;

	Model * model = new InvertedPendulum();
	InvertedPendulum * modelspecify = ( ( InvertedPendulum * ) model );
	ModelShootingData * shooting = new InvertedPendulumShootingData ( modelspecify );
	Simulator * simulator = new Simulator ( model );

	OdeManager * odemanager = new CacheOdeManager();
	yane::MinProg::NLP * sqp = new SqpFortran();

	SqpFortran * sqpSpecify = ( ( SqpFortran * ) sqp );
	sqpSpecify->setAccuracy ( 1E-4 );
	sqpSpecify->setMaxFun ( 20 );
	sqpSpecify->setMaxIterations ( 1000 );
	sqpSpecify->setLineSearchTol ( 0.1 );

	MPC * mpc = new MPC ( INF );
	mpc->reset ( odemanager, sqp, model, HORIZONMAX, PRECEDESTEPS, shooting );

	double computingtime, closedloopcost;
	double * nextControls = new double [ model->dimensionControl() * HORIZONSHIFT ];
	double * nextTimesteps = new double [ HORIZONSHIFT + 1 ];
	double * valuefunction = new double [ HORIZONMAX ];
	double * trajectory = new double [ ( HORIZONMAX + 1 ) * model->dimensionState() ];
	double * odetolerance = new double [ 9 ];
	double * sqptolerance = new double [ 9 ];
	for ( int i = 0; i < 9; i++ )
	{
		odetolerance[i] = pow(0.1, i+1);
		sqptolerance[i] = pow(0.1, i+1);
	}

	for ( int k = 0; k < 9; k++ )
	{
		for ( int l = 0; l < 9; l++ )
		{
			sqpSpecify->setAccuracy ( sqptolerance[k] );
			mpc->allocateMemory ( t, u );
			mpc->generateTimeGrid ( t, 0.0, H_NEW );
			mpc->initCalc ( t, u );
			mpc->resizeHorizon ( optimizationHorizon, H_NEW );
			mpc->setConfigShiftHorizon ( STANDARD, CURRENT_DATA );
			x = new double [ model->dimensionState() ];
			model->getDefaultState ( x );

			closedloopcost = 0.0;
			cout << "Optimization tolerance: " << setw(7) << setprecision ( 5 ) << sqptolerance[k] << "  Ode solver tolerance: " << setw(7) << setprecision ( 5 ) << odetolerance[l] << endl;
			cout << "  Time          u        phi      omega          x          v      ms        Cost" << endl;
			cout << "=================================================================================" << endl;
			for ( int i = 0; i < 200; i++ )
			{
				RTClock timer;
				try
				{
					model->odeConfig()->setTolerance ( odetolerance[l], odetolerance[l] );
					mpc->calc ( x );
				}
				catch ( yane::MinProg::SolverWarning e )
				{
				}
				double duration = timer.elapsedSeconds() * 1000.0;

				odemanager->getObjFunValues ( t, x, u, valuefunction, trajectory );
				closedloopcost += valuefunction[0];

				try
				{
					mpc->shiftHorizon ( nextControls, nextTimesteps, HORIZONSHIFT, x, H_NEW );
				}
				catch ( yane::MinProg::SolverWarning e )
				{
				}
				for ( int j = 0; j < HORIZONSHIFT; j++ )
				{
					cout << fixed << setprecision ( 3 ) << setw ( 6 ) << nextTimesteps[j] << " " << setprecision ( 5 ) << setw ( 10 ) << nextControls[j] << " ";
					for ( int i = 0; i < 4; i++ )
					{
						cout << setprecision ( 5 ) << setw ( 10 ) << x[i] << " ";
					}
					if ( j == 0 )
					{
						cout << setw ( 7 ) << setprecision ( 2 ) << duration;
					}
					else
					{
						cout << "     --";
					}
					cout << setw ( 12 ) << setprecision ( 3 ) << closedloopcost;
					cout << endl;

					modelspecify->odeConfig()->setTolerance ( 1e-10, 1e-10 );
					simulator->initialize ( nextTimesteps[j], x, &nextControls[j] );
					simulator->run ( nextTimesteps[j+1] );
				}
			}
			delete[] x;
		}
	}

	delete mpc;
	delete odemanager;
	delete sqp;
	delete simulator;
	delete shooting;
	delete model;
	delete[] t;
	delete[] u;
	delete[] nextControls;
	delete[] nextTimesteps;
	delete[] trajectory;
	delete[] valuefunction;
	delete[] odetolerance;
	delete[] sqptolerance;

	return 0;
}
