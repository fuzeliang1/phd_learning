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

#include "model_7_14.h"

#include <yane/odesolve.h>
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

#define PI 3.141592654
#define INF 1.0E19

const yane::Utils::Uuid yane::ExampleModels2::InvertedPendulum2::ID (
        ( const unsigned char* ) ( "507df855-d1fd-4470-831d-92eaef4f71c1" ) );

yane::ExampleModels2::InvertedPendulum2::InvertedPendulum2 ( ) :
	yane::Model::Model ( new yane::OdeSolve::DoPri853 ( ),
	        new yane::OdeSolve::DoPriConfig ( ), 4, 1, 2, 3 )
{
	setDoubleParameter ( 0, 0.01 );
	setDoubleParameter ( 1, 10.0 );
	setDoubleParameter ( 2, 0.01 );

	odeConfig ( )->setTolerance ( 1E-7, 1E-7 );

	_state_lb [ 0 ] = - 1E19;
	_state_lb [ 1 ] = - 1E19;
	_state_lb [ 2 ] = - 1E19;
	_state_lb [ 3 ] = - 1E19;

	_state_ub [ 0 ] = 1E19;
	_state_ub [ 1 ] = 1E19;
	_state_ub [ 2 ] = 1E19;
	_state_ub [ 3 ] = 1E19;

	_control_lb [ 0 ] = - 15.0;
	_control_ub [ 0 ] = 15.0;

	_terminalTime = new double [ 1 ];
}

yane::ExampleModels2::InvertedPendulum2::~InvertedPendulum2 ( )
{
	delete odeSolver ( );
	delete odeConfig ( );
	delete [ ] _terminalTime;
}

void yane::ExampleModels2::InvertedPendulum2::odeFunction ( double t,
        double * x, double * u, double * dx )
{
	double reibung = ( 4.0 * x [ 1 ] / ( 1.0 + 4.0 * x [ 1 ] * x [ 1 ] ) + 2.0
	        * atan ( 2.0 * x [ 1 ] ) / PI ) * _params [ 2 ];

	dx [ 0 ] = x [ 1 ];
	dx [ 1 ] = - ( 9.81 * sin ( x [ 0 ] ) + u [ 0 ] * cos ( x [ 0 ] ) )
	        / _params [ 1 ] - _params [ 0 ] * x [ 1 ] * abs ( x [ 1 ] )
	        - reibung;
	dx [ 2 ] = x [ 3 ];
	dx [ 3 ] = u [ 0 ];
}

double yane::ExampleModels2::InvertedPendulum2::continuousCostfunction (
        double t, double * x, double * u )
{
	double sinxpi = sin ( x [ 0 ] - PI );
	double cosy = cos ( x [ 1 ] );
	double temp = ( 1.0 - cos ( x [ 0 ] - PI ) ) * ( 1.0 + cosy * cosy );

	if ( t < _terminalTime [ 0 ] )
	{
		return pow ( 3.51 * sinxpi * sinxpi + ( 4.82 * sinxpi + 2.31 * x [ 1 ] )
		        * x [ 1 ] + 0.1 * temp * temp + 0.01 * x [ 2 ] * x [ 2 ] + 0.1
		        * x [ 3 ] * x [ 3 ], 2.0 ) + 1.0e-4 * u [ 0 ] * u [ 0 ];
	}
	else
	{
		return _terminalWeight * ( pow ( 3.51 * sinxpi * sinxpi + ( 4.82
		        * sinxpi + 2.31 * x [ 1 ] ) * x [ 1 ] + 0.1 * temp * temp
		        + 0.01 * x [ 2 ] * x [ 2 ] + 0.1 * x [ 3 ] * x [ 3 ], 2.0 )
		        + 1.0e-4 * u [ 0 ] * u [ 0 ] );
	}
}

void yane::ExampleModels2::InvertedPendulum2::restrictionFunction ( double t,
        double * x, double * u, double * fx )
{
	fx [ 0 ] = - x [ 0 ] + PI;
	fx [ 1 ] = x [ 0 ] + 3.0 * PI;
}

double yane::ExampleModels2::InvertedPendulum2::discreteCostfunction (
        int length, int horizont, double * t, double * x, double * u )
{
}

void yane::ExampleModels2::InvertedPendulum2::getObjectiveWeight (
        double & weight_continuousCostfunction,
        double & weight_discreteCostfunction )
{
	weight_continuousCostfunction = 1.0;
	weight_discreteCostfunction = 0.0;
}

void yane::ExampleModels2::InvertedPendulum2::getDefaultState ( double * x )
{
	x [ 0 ] = 1.5;
	x [ 1 ] = 0.0;
	x [ 2 ] = 0.0;
	x [ 3 ] = 0.0;
}

void yane::ExampleModels2::InvertedPendulum2::getDefaultControl ( double * u )
{
	u [ 0 ] = - 10.0;
}

void yane::ExampleModels2::InvertedPendulum2::eventBeforeMPC ( int horizont,
        double * t, double * x, double * sdatavalues )
{
}

void yane::ExampleModels2::InvertedPendulum2::setTerminalWeight (
        double terminalWeight, double * terminalTime )
{
	_terminalWeight = terminalWeight;
	_terminalTime = terminalTime;
}

double yane::ExampleModels2::InvertedPendulum2::getTerminalWeight ( void )
{
	return _terminalWeight;
}

const yane::Utils::Uuid & yane::ExampleModels2::InvertedPendulum2::id ( )
{
	return ID;
}
