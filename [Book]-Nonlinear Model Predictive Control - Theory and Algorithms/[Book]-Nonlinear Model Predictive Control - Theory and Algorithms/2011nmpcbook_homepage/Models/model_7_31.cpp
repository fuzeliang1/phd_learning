/***************************************************************************
 *                                                                         *
 * Universitaet Bayreuth - Lehrstuhl fuer Angewandte Mathematik            *
 *                                                                         *
 * Authors:                                                                *
 *   Juergen Pannek <juergen.pannek@uni-bayreuth.de>                       *
 *                                                                         *
 ***************************************************************************/

/**
 @file
 @author Juergen Pannek <juergen.pannek@uni-bayreuth.de>
 */

#include "model_7_31.h"

#include <yane/odesolve.h>
#include <iostream>
#include <cmath>

using namespace std;

#define PI 3.141592654
#define INF 1.0E19

const yane::Utils::Uuid yane::ExampleModels2::Crane::ID (
        ( const unsigned char* ) ( "efed7a0b-f8d2-4e17-993e-7dc9e2642b61" ) );

yane::ExampleModels2::Crane::Crane ( ) :
	yane::Model::Model ( new yane::OdeSolve::DoPri853 ( ),
	        new yane::OdeSolve::DoPriConfig ( ), 4, 1, 0, 4 )
{
	setDoubleParameter ( 0, 9.81 );
	setDoubleParameter ( 1, 1.25 );
	setDoubleParameter ( 2, 0.007 );
	setDoubleParameter ( 3, 0.197 );

	odeConfig ( )->setTolerance ( 1E-8, 1E-8 );

	_state_lb [ 0 ] = - INF;
	_state_lb [ 1 ] = - INF;
	_state_lb [ 2 ] = - INF;
	_state_lb [ 3 ] = - INF;

	_state_ub [ 0 ] = INF;
	_state_ub [ 1 ] = INF;
	_state_ub [ 2 ] = INF;
	_state_ub [ 3 ] = INF;

	_control_lb [ 0 ] = - 1.0;
	_control_ub [ 0 ] = 1.0;

}

yane::ExampleModels2::Crane::~Crane ( )
{
	delete odeSolver ( );
	delete odeConfig ( );
}

void yane::ExampleModels2::Crane::odeFunction ( double t, double * x,
        double * u, double * dx )
{
	double friction = ( 4.0 * x [ 1 ] / ( 1.0 + 4.0 * x [ 1 ] * x [ 1 ] ) + 2.0
	        * atan ( 2.0 * x [ 1 ] ) / PI ) * _params [ 3 ];
	double cubic = _params [ 2 ] / _params [ 1 ] * x [ 1 ] * x [ 1 ] * atan (
	        1.0e3 * x [ 1 ] );

	dx [ 0 ] = x [ 1 ];
	dx [ 1 ] = - _params [ 0 ] / _params [ 1 ] * sin ( x [ 0 ] ) - u [ 0 ]
	        * cos ( x [ 0 ] ) - cubic - friction;
	dx [ 2 ] = x [ 3 ];
	dx [ 3 ] = u [ 0 ];
}

double yane::ExampleModels2::Crane::continuousCostfunction ( double t,
        double * x, double * u )
{
	return 100.0 * sin ( 0.5 * x [ 0 ] ) * sin ( 0.5 * x [ 0 ] ) + x [ 1 ]
	        * x [ 1 ] + 10.0 * x [ 2 ] * x [ 2 ] + x [ 3 ] * x [ 3 ] + u [ 0 ]
	        * u [ 0 ];
}

void yane::ExampleModels2::Crane::restrictionFunction ( double t, double * x,
        double * u, double * fx )
{
}

double yane::ExampleModels2::Crane::discreteCostfunction ( int length,
        int horizon, double * t, double * x, double * u )
{
	return 0.0;
}

void yane::ExampleModels2::Crane::getObjectiveWeight (
        double & weight_continuousCostfunction,
        double & weight_discreteCostfunction )
{
	weight_continuousCostfunction = 1.0;
	weight_discreteCostfunction = 0.0;
}

void yane::ExampleModels2::Crane::getDefaultState ( double * x )
{
	x [ 0 ] = 0.0;
	x [ 1 ] = 0.0;
	x [ 2 ] = 10.0;
	x [ 3 ] = 0.0;
}

void yane::ExampleModels2::Crane::getDefaultControl ( double * u )
{
	u [ 0 ] = 0.0;
}

void yane::ExampleModels2::Crane::eventBeforeMPC ( int horizon,
        double * timesteps, double * x, double * sdatavalues )
{
}

const yane::Utils::Uuid & yane::ExampleModels2::Crane::id ( )
{
	return ID;
}
