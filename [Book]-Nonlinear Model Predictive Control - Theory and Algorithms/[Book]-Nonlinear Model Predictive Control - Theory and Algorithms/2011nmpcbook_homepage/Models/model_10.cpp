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

#include "model_10.h"

#include <yane/odesolve.h>
#include <iostream>
#include <cmath>

using namespace std;

#define PI 3.141592654
#define INF 1.0E19

const yane::Utils::Uuid yane::ExampleModels2::InvertedPendulum::ID (
        ( const unsigned char* ) ( "465333aa-245d-47fd-86d1-a6eb73d671f0" ) );

yane::ExampleModels2::InvertedPendulum::InvertedPendulum ( ) :
	yane::Model::Model ( new yane::OdeSolve::DoPri853 ( ),
	        new yane::OdeSolve::DoPriConfig ( ), 4, 1, 0, 4 )
{
	setDoubleParameter ( 0, 9.81 );
	setDoubleParameter ( 1, 1.25 );
	setDoubleParameter ( 2, 0.007 );
	setDoubleParameter ( 3, 0.197 );

	odeConfig ( )->setTolerance ( 1E-10, 1E-10 );

	_state_lb [ 0 ] = - 1E19;
	_state_lb [ 1 ] = - 1E19;
	_state_lb [ 2 ] = - 5.0;
	_state_lb [ 3 ] = - 10.0;

	_state_ub [ 0 ] = 1E19;
	_state_ub [ 1 ] = 1E19;
	_state_ub [ 2 ] = 5.0;
	_state_ub [ 3 ] = 10.0;

	_control_lb [ 0 ] = - 5.0;
	_control_ub [ 0 ] = 5.0;
}

yane::ExampleModels2::InvertedPendulum::~InvertedPendulum ( )
{
	delete odeSolver ( );
	delete odeConfig ( );
}

void yane::ExampleModels2::InvertedPendulum::odeFunction ( double t,
        double * x, double * u, double * dx )
{
	dx [ 0 ] = x [ 1 ];
	dx [ 1 ] = - _params [ 0 ] / _params [ 1 ] * sin ( x [ 0 ] ) - _params [ 2 ]
	        / _params [ 1 ] * atan ( 1e3 * x [ 1 ] ) * x [ 1 ] * x [ 1 ]
	        - u [ 0 ] * cos ( x [ 0 ] ) - _params [ 3 ] * ( 4.0 * x [ 1 ]
	        / ( 1.0 + 4.0 * x [ 1 ] * x [ 1 ] ) + 2.0 * atan ( 2.0 * x [ 1 ] )
	        / PI );
	dx [ 2 ] = x [ 3 ];
	dx [ 3 ] = u [ 0 ];
}

double yane::ExampleModels2::InvertedPendulum::continuousCostfunction (
        double t, double * x, double * u )
{
	double sinxpi = sin ( x [ 0 ] - PI );
	double cosy = cos ( x [ 1 ] );
	double temp = ( 1.0 - cos ( x [ 0 ] - PI ) ) * ( 1.0 + cosy * cosy );

	return pow ( 3.51 * sinxpi * sinxpi + ( 4.82 * sinxpi + 2.31 * x [ 1 ] )
	        * x [ 1 ] + 2.0 * temp * temp + x [ 2 ] * x [ 2 ] + x [ 3 ]
	        * x [ 3 ] + u [ 0 ] * u [ 0 ], 2.0 );
}

void yane::ExampleModels2::InvertedPendulum::restrictionFunction ( double t,
        double * x, double * u, double * fx )
{
}

double yane::ExampleModels2::InvertedPendulum::discreteCostfunction (
        int length, int horizon, double * t, double * x, double * u )
{
	return 0.0;
}

void yane::ExampleModels2::InvertedPendulum::getObjectiveWeight (
        double & weight_continuousCostfunction,
        double & weight_discreteCostfunction )
{
	weight_continuousCostfunction = 1.0;
	weight_discreteCostfunction = 0.0;
}

void yane::ExampleModels2::InvertedPendulum::getDefaultState ( double * x )
{
	x [ 0 ] = 2.0;
	x [ 1 ] = 2.0;
	x [ 2 ] = 0.0;
	x [ 3 ] = 0.0;
}

void yane::ExampleModels2::InvertedPendulum::getDefaultControl ( double * u )
{
	u [ 0 ] = 0.0;
}

void yane::ExampleModels2::InvertedPendulum::eventBeforeMPC ( int horizon,
        double * timesteps, double * x, double * sdatavalues )
{
}

const yane::Utils::Uuid & yane::ExampleModels2::InvertedPendulum::id ( )
{
	return ID;
}
