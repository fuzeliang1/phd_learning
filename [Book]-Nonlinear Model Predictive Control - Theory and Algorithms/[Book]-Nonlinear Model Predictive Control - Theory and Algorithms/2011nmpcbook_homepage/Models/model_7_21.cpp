/***************************************************************************
 *                                                                         *
 * Arm-Rotor-Platform Example                                              *
 *                                                                         *
 * Authors:                                                                *
 *   Juergen Pannek <juergen.pannek@uni-bayreuth.de>                       *
 *                                                                         *
 ***************************************************************************/

/**
 @file
 @author Juergen Pannek <juergen.pannek@uni-bayreuth.de>
 */

#include "model_7_21.h"

#include <yane/odesolve.h>
#include <iostream>
#include <cmath>

using namespace std;

#define PI 3.141592654
#define INF 1.0E19

const yane::Utils::Uuid yane::ExampleModels2::ARP::ID (
        ( const unsigned char* ) ( "710eb51f-3364-4d5c-9eef-12a6e7233f38" ) );

yane::ExampleModels2::ARP::ARP ( ) :
	yane::Model::Model ( new yane::OdeSolve::DoPri853 ( ),
	        new yane::OdeSolve::DoPriConfig ( ), 16, 1, 0, 28 )
{
	odeConfig ( )->setTolerance ( 1E-10, 1E-10 );

	setModelParameter ( );

	_state_lb [ 0 ] = - INF;
	_state_lb [ 1 ] = - INF;
	_state_lb [ 2 ] = - INF;
	_state_lb [ 3 ] = - INF;
	_state_lb [ 4 ] = - INF;
	_state_lb [ 5 ] = - INF;
	_state_lb [ 6 ] = - INF;
	_state_lb [ 7 ] = - INF;
	_state_lb [ 8 ] = - INF;
	_state_lb [ 9 ] = - INF;
	_state_lb [ 10 ] = - INF;
	_state_lb [ 11 ] = - INF;
	_state_lb [ 12 ] = - INF;
	_state_lb [ 13 ] = - INF;
	_state_lb [ 14 ] = - INF;
	_state_lb [ 15 ] = - INF;

	_state_ub [ 0 ] = INF;
	_state_ub [ 1 ] = INF;
	_state_ub [ 2 ] = INF;
	_state_ub [ 3 ] = INF;
	_state_ub [ 4 ] = INF;
	_state_ub [ 5 ] = INF;
	_state_ub [ 6 ] = INF;
	_state_ub [ 7 ] = INF;
	_state_ub [ 8 ] = INF;
	_state_ub [ 9 ] = INF;
	_state_ub [ 10 ] = INF;
	_state_ub [ 11 ] = INF;
	_state_ub [ 12 ] = INF;
	_state_ub [ 13 ] = INF;
	_state_ub [ 14 ] = INF;
	_state_ub [ 15 ] = INF;

	_control_lb [ 0 ] = - INF;
	_control_ub [ 0 ] = INF;
}

yane::ExampleModels2::ARP::~ARP ( )
{
	delete odeSolver ( );
	delete odeConfig ( );
}

void yane::ExampleModels2::ARP::odeFunction ( double t, double * x, double * u,
        double * dx )
{
	double zetaref0 = _params [ 22 ]; // Reference Function
	double zetaref1 = _params [ 23 ];
	double zetaref2 = _params [ 24 ];
	double zetaref3 = _params [ 25 ];
	double zetaref4 = _params [ 26 ];

	double
	        v =
	                - _params [ 18 ] * ( x [ 8 ]
	                        - _params [ 12 ] / ( _params [ 10 ] - _params [ 11 ]
	                                * _params [ 12 ] ) * ( x [ 10 ]
	                                - _params [ 12 ] * x [ 12 ] ) - zetaref0 )
	                        - _params [ 19 ] * ( _params [ 10 ]
	                                / ( _params [ 10 ] - _params [ 11 ]
	                                        * _params [ 12 ] ) * x [ 10 ]
	                                - _params [ 12 ] / ( _params [ 10 ]
	                                        - _params [ 11 ] * _params [ 12 ] )
	                                        * ( - _params [ 10 ] * x [ 8 ]
	                                                + _params [ 10 ] * x [ 12 ]
	                                                - _params [ 16 ] * x [ 0 ]
	                                                - _params [ 17 ] * x [ 2 ] )
	                                - zetaref1 )
	                        - _params [ 20 ]
	                                * ( _params [ 10 ] / ( _params [ 10 ]
	                                        - _params [ 11 ] * _params [ 12 ] )
	                                        * ( - _params [ 10 ] * x [ 8 ]
	                                                - _params [ 11 ] * x [ 10 ]
	                                                + _params [ 10 ] * x [ 12 ]
	                                                - _params [ 16 ] * x [ 0 ]
	                                                - _params [ 17 ] * x [ 2 ] )
	                                        - _params [ 12 ] / ( _params [ 10 ]
	                                                - _params [ 11 ]
	                                                        * _params [ 12 ] )
	                                                * ( - _params [ 10 ]
	                                                        * x [ 10 ]
	                                                        - _params [ 16 ]
	                                                                * ( x [ 2 ]
	                                                                        + x [ 4 ]
	                                                                                * x [ 10 ] )
	                                                        - _params [ 17 ]
	                                                                * ( - _params [ 5 ]
	                                                                        / _params [ 0 ]
	                                                                        * x [ 0 ]
	                                                                        - _params [ 7 ]
	                                                                                / _params [ 0 ]
	                                                                                * x [ 2 ]
	                                                                        + x [ 10 ]
	                                                                                * x [ 6 ]
	                                                                        - _params [ 7 ]
	                                                                                * _params [ 1 ]
	                                                                                * _params [ 2 ]
	                                                                                / pow (
	                                                                                        _params [ 0 ],
	                                                                                        2.0 )
	                                                                                * x [ 10 ] ) )
	                                        - zetaref2 )
	                        - _params [ 21 ]
	                                * ( 1.0 / ( _params [ 10 ] - _params [ 11 ]
	                                        * _params [ 12 ] )
	                                        * ( - pow ( _params [ 10 ], 2.0 )
	                                                * ( x [ 10 ] - x [ 14 ] )
	                                                + ( - _params [ 10 ]
	                                                        * x [ 8 ]
	                                                        - _params [ 11 ]
	                                                                * x [ 10 ]
	                                                        + _params [ 10 ]
	                                                                * x [ 12 ]
	                                                        + _params [ 12 ]
	                                                                * x [ 14 ]
	                                                        - _params [ 16 ]
	                                                                * x [ 0 ]
	                                                        - _params [ 17 ]
	                                                                * x [ 2 ] )
	                                                        * ( - _params [ 10 ]
	                                                                * ( _params [ 11 ]
	                                                                        - _params [ 12 ] )
	                                                                + _params [ 12 ]
	                                                                        * ( _params [ 16 ]
	                                                                                * x [ 4 ]
	                                                                                + _params [ 17 ]
	                                                                                        * ( x [ 6 ]
	                                                                                                - _params [ 7 ]
	                                                                                                        * _params [ 1 ]
	                                                                                                        * _params [ 2 ]
	                                                                                                        / pow (
	                                                                                                                _params [ 0 ],
	                                                                                                                2.0 ) ) ) )
	                                                + ( ( - _params [ 10 ]
	                                                        * _params [ 16 ]
	                                                        - _params [ 12 ]
	                                                                * _params [ 17 ]
	                                                                * _params [ 5 ]
	                                                                / _params [ 0 ] )
	                                                        * ( - _params [ 5 ]
	                                                                / _params [ 0 ]
	                                                                * x [ 0 ]
	                                                                - _params [ 7 ]
	                                                                        / _params [ 0 ]
	                                                                        * x [ 2 ]
	                                                                + x [ 10 ]
	                                                                        * x [ 6 ]
	                                                                - _params [ 7 ]
	                                                                        * _params [ 1 ]
	                                                                        * _params [ 2 ]
	                                                                        / pow (
	                                                                                _params [ 0 ],
	                                                                                2.0 )
	                                                                        * x [ 10 ]
	                                                                + x [ 10 ]
	                                                                        * ( - x [ 0 ]
	                                                                                * x [ 10 ]
	                                                                                + x [ 6 ] ) )
	                                                        + ( - _params [ 10 ]
	                                                                * _params [ 17 ]
	                                                                + _params [ 12 ]
	                                                                        * _params [ 16 ]
	                                                                - _params [ 12 ]
	                                                                        * _params [ 7 ]
	                                                                        * _params [ 17 ]
	                                                                        / _params [ 0 ] )
	                                                                * ( - _params [ 5 ]
	                                                                        / _params [ 0 ]
	                                                                        * ( x [ 2 ]
	                                                                                + x [ 4 ]
	                                                                                        * x [ 10 ] )
	                                                                        - _params [ 7 ]
	                                                                                / _params [ 0 ]
	                                                                                * ( - _params [ 5 ]
	                                                                                        / _params [ 0 ]
	                                                                                        * x [ 0 ]
	                                                                                        - _params [ 7 ]
	                                                                                                / _params [ 0 ]
	                                                                                                * x [ 2 ]
	                                                                                        + x [ 10 ]
	                                                                                                * x [ 6 ]
	                                                                                        - _params [ 7 ]
	                                                                                                * _params [ 1 ]
	                                                                                                * _params [ 2 ]
	                                                                                                * x [ 10 ]
	                                                                                                / pow (
	                                                                                                        _params [ 0 ],
	                                                                                                        2.0 ) )
	                                                                        + x [ 10 ]
	                                                                                * ( - x [ 10 ]
	                                                                                        * x [ 2 ]
	                                                                                        - _params [ 5 ]
	                                                                                                / _params [ 0 ]
	                                                                                                * x [ 4 ]
	                                                                                        - _params [ 7 ]
	                                                                                                / _params [ 0 ]
	                                                                                                * x [ 6 ]
	                                                                                        + _params [ 5 ]
	                                                                                                * _params [ 1 ]
	                                                                                                * _params [ 2 ]
	                                                                                                / pow (
	                                                                                                        _params [ 0 ],
	                                                                                                        2.0 ) ) )
	                                                        + ( _params [ 12 ]
	                                                                * _params [ 16 ]
	                                                                * x [ 10 ] )
	                                                                * ( - x [ 10 ]
	                                                                        * ( x [ 2 ]
	                                                                                + x [ 10 ]
	                                                                                        * x [ 4 ] )
	                                                                        - x [ 2 ]
	                                                                                * x [ 10 ]
	                                                                        - _params [ 5 ]
	                                                                                / _params [ 0 ]
	                                                                                * x [ 4 ]
	                                                                        - _params [ 7 ]
	                                                                                / _params [ 0 ]
	                                                                                * x [ 6 ]
	                                                                        + _params [ 5 ]
	                                                                                * _params [ 1 ]
	                                                                                * _params [ 2 ]
	                                                                                / pow (
	                                                                                        _params [ 0 ],
	                                                                                        2.0 ) )
	                                                        + ( _params [ 12 ]
	                                                                * _params [ 17 ]
	                                                                * x [ 10 ] )
	                                                                * ( - x [ 10 ]
	                                                                        * ( - _params [ 5 ]
	                                                                                / _params [ 0 ]
	                                                                                * x [ 0 ]
	                                                                                - _params [ 7 ]
	                                                                                        / _params [ 0 ]
	                                                                                        * x [ 2 ]
	                                                                                + x [ 10 ]
	                                                                                        * x [ 6 ] )
	                                                                        - _params [ 5 ]
	                                                                                / _params [ 0 ]
	                                                                                * ( - x [ 0 ]
	                                                                                        * x [ 10 ]
	                                                                                        + x [ 6 ] )
	                                                                        - _params [ 7 ]
	                                                                                / _params [ 0 ]
	                                                                                * ( - x [ 10 ]
	                                                                                        * x [ 2 ]
	                                                                                        - _params [ 5 ]
	                                                                                                / _params [ 0 ]
	                                                                                                * x [ 4 ]
	                                                                                        - _params [ 7 ]
	                                                                                                / _params [ 0 ]
	                                                                                                * x [ 6 ]
	                                                                                        + _params [ 5 ]
	                                                                                                * _params [ 1 ]
	                                                                                                * _params [ 2 ]
	                                                                                                / pow (
	                                                                                                        _params [ 0 ],
	                                                                                                        2.0 ) ) ) ) )
	                                        - zetaref3 ) + zetaref4;

	dx [ 0 ] = x [ 2 ] + x [ 10 ] * x [ 4 ];
	dx [ 1 ] = x [ 3 ] + x [ 11 ] * x [ 5 ];

	dx [ 2 ] = - _params [ 5 ] / _params [ 0 ] * x [ 0 ] - _params [ 7 ]
	        / _params [ 0 ] * x [ 2 ] + x [ 10 ] * x [ 6 ] - ( _params [ 1 ]
	        * _params [ 2 ] ) / pow ( _params [ 0 ], 2.0 ) * _params [ 7 ]
	        * x [ 10 ];
	dx [ 3 ] = - _params [ 5 ] / _params [ 0 ] * x [ 1 ] - _params [ 7 ]
	        / _params [ 0 ] * x [ 3 ] + x [ 11 ] * x [ 7 ] - ( _params [ 1 ]
	        * _params [ 2 ] ) / pow ( _params [ 0 ], 2.0 ) * _params [ 7 ]
	        * x [ 11 ];

	dx [ 4 ] = - x [ 10 ] * x [ 0 ] + x [ 6 ];
	dx [ 5 ] = - x [ 11 ] * x [ 1 ] + x [ 7 ];

	dx [ 6 ] = - x [ 10 ] * x [ 2 ] - _params [ 5 ] / _params [ 0 ] * x [ 4 ]
	        - _params [ 7 ] / _params [ 0 ] * x [ 6 ] + ( _params [ 1 ]
	        * _params [ 2 ] ) / pow ( _params [ 0 ], 2.0 ) * _params [ 5 ];
	dx [ 7 ] = - x [ 11 ] * x [ 3 ] - _params [ 5 ] / _params [ 0 ] * x [ 5 ]
	        - _params [ 7 ] / _params [ 0 ] * x [ 7 ] + ( _params [ 1 ]
	        * _params [ 2 ] ) / pow ( _params [ 0 ], 2.0 ) * _params [ 5 ];

	dx [ 8 ] = x [ 10 ];
	dx [ 9 ] = x [ 11 ];

	dx [ 10 ] = - _params [ 10 ] * x [ 8 ] - _params [ 11 ] * x [ 10 ]
	        + _params [ 10 ] * x [ 12 ] + _params [ 12 ] * x [ 14 ]
	        - _params [ 16 ] * x [ 0 ] - _params [ 17 ] * x [ 2 ];
	dx [ 11 ] = - _params [ 10 ] * x [ 9 ] - _params [ 11 ] * x [ 11 ]
	        + _params [ 10 ] * x [ 13 ] + _params [ 12 ] * x [ 15 ]
	        - _params [ 16 ] * x [ 1 ] - _params [ 17 ] * x [ 3 ];

	dx [ 12 ] = x [ 14 ];
	dx [ 13 ] = x [ 15 ];

	dx [ 14 ]
	        = _params [ 13 ] * x [ 8 ] + _params [ 14 ] * x [ 10 ]
	                - _params [ 13 ] * x [ 12 ] - ( _params [ 14 ]
	                + _params [ 15 ] ) * x [ 14 ]
	                + 1.0 / _params [ 4 ]
	                        * ( _params [ 4 ] / ( pow ( _params [ 10 ], 2.0 )
	                                + _params [ 12 ] * ( _params [ 16 ]
	                                        * x [ 4 ] + _params [ 17 ]
	                                        * ( x [ 6 ] - _params [ 7 ]
	                                                * _params [ 1 ]
	                                                * _params [ 2 ] / pow (
	                                                _params [ 0 ], 2.0 ) ) ) )
	                                * ( - ( - _params [ 10 ] * x [ 8 ]
	                                        - _params [ 11 ] * x [ 10 ]
	                                        + _params [ 10 ] * x [ 12 ]
	                                        + _params [ 12 ] * x [ 14 ]
	                                        - _params [ 16 ] * x [ 0 ]
	                                        - _params [ 17 ] * x [ 2 ] )
	                                        * ( - pow ( _params [ 10 ], 2.0 )
	                                                + _params [ 10 ]
	                                                        * _params [ 11 ]
	                                                        * ( _params [ 11 ]
	                                                                - _params [ 12 ] )
	                                                + ( - ( _params [ 10 ]
	                                                        + _params [ 11 ]
	                                                                * _params [ 12 ] )
	                                                        * _params [ 16 ]
	                                                        - _params [ 5 ]
	                                                                / _params [ 0 ]
	                                                                * _params [ 12 ]
	                                                                * _params [ 17 ] )
	                                                        * x [ 4 ]
	                                                + ( - ( _params [ 10 ]
	                                                        + _params [ 11 ]
	                                                                * _params [ 12 ] )
	                                                        * _params [ 17 ]
	                                                        + _params [ 12 ]
	                                                                * _params [ 16 ]
	                                                        - _params [ 7 ]
	                                                                / _params [ 0 ]
	                                                                * _params [ 12 ]
	                                                                * _params [ 17 ] )
	                                                        * ( x [ 6 ]
	                                                                - _params [ 7 ]
	                                                                        * _params [ 1 ]
	                                                                        * _params [ 2 ]
	                                                                        / pow (
	                                                                                _params [ 0 ],
	                                                                                2.0 ) )
	                                                + _params [ 12 ]
	                                                        * _params [ 16 ]
	                                                        * x [ 10 ]
	                                                        * ( - x [ 0 ] )
	                                                + _params [ 12 ]
	                                                        * _params [ 17 ]
	                                                        * x [ 10 ]
	                                                        * ( - x [ 2 ] )
	                                                + 2.0 * _params [ 12 ]
	                                                        * _params [ 16 ]
	                                                        * ( - x [ 10 ]
	                                                                * x [ 0 ]
	                                                                + x [ 6 ] )
	                                                + 2.0 * _params [ 12 ]
	                                                        * _params [ 17 ]
	                                                        * ( - x [ 10 ]
	                                                                * x [ 2 ]
	                                                                - _params [ 5 ]
	                                                                        / _params [ 0 ]
	                                                                        * x [ 4 ]
	                                                                - _params [ 7 ]
	                                                                        / _params [ 0 ]
	                                                                        * x [ 6 ] ) )
	                                        - ( _params [ 13 ] * x [ 8 ]
	                                                + _params [ 14 ] * x [ 10 ]
	                                                - _params [ 13 ] * x [ 12 ]
	                                                - ( _params [ 14 ]
	                                                        + _params [ 15 ] )
	                                                        * x [ 14 ] )
	                                                * ( pow ( _params [ 10 ],
	                                                        2.0 )
	                                                        + _params [ 12 ]
	                                                                * ( - _params [ 10 ]
	                                                                        * ( _params [ 11 ]
	                                                                                - _params [ 12 ] )
	                                                                        + _params [ 12 ]
	                                                                                * ( _params [ 16 ]
	                                                                                        * x [ 4 ]
	                                                                                        + _params [ 17 ]
	                                                                                                * ( x [ 6 ]
	                                                                                                        - _params [ 7 ]
	                                                                                                                * _params [ 1 ]
	                                                                                                                * _params [ 2 ]
	                                                                                                                / pow (
	                                                                                                                        _params [ 0 ],
	                                                                                                                        2.0 ) ) ) ) )
	                                        - ( ( - _params [ 10 ]
	                                                * _params [ 16 ]
	                                                - _params [ 12 ]
	                                                        * _params [ 17 ]
	                                                        * _params [ 5 ]
	                                                        / _params [ 0 ] )
	                                                * ( - _params [ 5 ]
	                                                        / _params [ 0 ]
	                                                        * x [ 0 ]
	                                                        - _params [ 7 ]
	                                                                / _params [ 0 ]
	                                                                * x [ 2 ]
	                                                        + x [ 10 ] * x [ 6 ]
	                                                        - _params [ 7 ]
	                                                                * _params [ 1 ]
	                                                                * _params [ 2 ]
	                                                                / pow (
	                                                                        _params [ 0 ],
	                                                                        2.0 )
	                                                                * x [ 10 ]
	                                                        + x [ 10 ]
	                                                                * ( - x [ 0 ]
	                                                                        * x [ 10 ]
	                                                                        + x [ 6 ] ) )
	                                                + ( - _params [ 10 ]
	                                                        * _params [ 17 ]
	                                                        + _params [ 12 ]
	                                                                * _params [ 16 ]
	                                                        - _params [ 12 ]
	                                                                * _params [ 7 ]
	                                                                * _params [ 17 ]
	                                                                / _params [ 0 ] )
	                                                        * ( - _params [ 5 ]
	                                                                / _params [ 0 ]
	                                                                * ( x [ 2 ]
	                                                                        + x [ 4 ]
	                                                                                * x [ 10 ] )
	                                                                - _params [ 7 ]
	                                                                        / _params [ 0 ]
	                                                                        * ( - _params [ 5 ]
	                                                                                / _params [ 0 ]
	                                                                                * x [ 0 ]
	                                                                                - _params [ 7 ]
	                                                                                        / _params [ 0 ]
	                                                                                        * x [ 2 ]
	                                                                                + x [ 10 ]
	                                                                                        * x [ 6 ]
	                                                                                - _params [ 7 ]
	                                                                                        * _params [ 1 ]
	                                                                                        * _params [ 2 ]
	                                                                                        * x [ 10 ]
	                                                                                        / pow (
	                                                                                                _params [ 0 ],
	                                                                                                2.0 ) )
	                                                                + x [ 10 ]
	                                                                        * ( - x [ 10 ]
	                                                                                * x [ 2 ]
	                                                                                - _params [ 5 ]
	                                                                                        / _params [ 0 ]
	                                                                                        * x [ 4 ]
	                                                                                - _params [ 7 ]
	                                                                                        / _params [ 0 ]
	                                                                                        * x [ 6 ]
	                                                                                + _params [ 5 ]
	                                                                                        * _params [ 1 ]
	                                                                                        * _params [ 2 ]
	                                                                                        / pow (
	                                                                                                _params [ 0 ],
	                                                                                                2.0 ) ) )
	                                                + ( _params [ 12 ]
	                                                        * _params [ 16 ]
	                                                        * x [ 10 ] )
	                                                        * ( - x [ 10 ]
	                                                                * ( x [ 2 ]
	                                                                        + x [ 10 ]
	                                                                                * x [ 4 ] )
	                                                                - x [ 2 ]
	                                                                        * x [ 10 ]
	                                                                - _params [ 5 ]
	                                                                        / _params [ 0 ]
	                                                                        * x [ 4 ]
	                                                                - _params [ 7 ]
	                                                                        / _params [ 0 ]
	                                                                        * x [ 6 ]
	                                                                + _params [ 5 ]
	                                                                        * _params [ 1 ]
	                                                                        * _params [ 2 ]
	                                                                        / pow (
	                                                                                _params [ 0 ],
	                                                                                2.0 ) )
	                                                + ( _params [ 12 ]
	                                                        * _params [ 17 ]
	                                                        * x [ 10 ] )
	                                                        * ( - x [ 10 ]
	                                                                * ( - _params [ 5 ]
	                                                                        / _params [ 0 ]
	                                                                        * x [ 0 ]
	                                                                        - _params [ 7 ]
	                                                                                / _params [ 0 ]
	                                                                                * x [ 2 ]
	                                                                        + x [ 10 ]
	                                                                                * x [ 6 ] )
	                                                                - _params [ 5 ]
	                                                                        / _params [ 0 ]
	                                                                        * ( - x [ 0 ]
	                                                                                * x [ 10 ]
	                                                                                + x [ 6 ] )
	                                                                - _params [ 7 ]
	                                                                        / _params [ 0 ]
	                                                                        * ( - x [ 10 ]
	                                                                                * x [ 2 ]
	                                                                                - _params [ 5 ]
	                                                                                        / _params [ 0 ]
	                                                                                        * x [ 4 ]
	                                                                                - _params [ 7 ]
	                                                                                        / _params [ 0 ]
	                                                                                        * x [ 6 ]
	                                                                                + _params [ 5 ]
	                                                                                        * _params [ 1 ]
	                                                                                        * _params [ 2 ]
	                                                                                        / pow (
	                                                                                                _params [ 0 ],
	                                                                                                2.0 ) ) ) )
	                                        - ( - _params [ 10 ] * ( x [ 10 ]
	                                                - x [ 14 ] )
	                                                - _params [ 16 ]
	                                                        * ( x [ 2 ]
	                                                                + x [ 4 ]
	                                                                        * x [ 10 ] )
	                                                - _params [ 17 ]
	                                                        * ( - _params [ 5 ]
	                                                                / _params [ 0 ]
	                                                                * x [ 0 ]
	                                                                - _params [ 7 ]
	                                                                        / _params [ 0 ]
	                                                                        * x [ 2 ]
	                                                                + x [ 6 ]
	                                                                        * x [ 10 ]
	                                                                - _params [ 7 ]
	                                                                        * _params [ 1 ]
	                                                                        * _params [ 2 ]
	                                                                        / pow (
	                                                                                _params [ 0 ],
	                                                                                2.0 )
	                                                                        * x [ 10 ] ) )
	                                                * ( - _params [ 10 ]
	                                                        * ( _params [ 11 ]
	                                                                - _params [ 12 ] )
	                                                        + _params [ 12 ]
	                                                                * _params [ 16 ]
	                                                                * x [ 4 ]
	                                                        + _params [ 17 ]
	                                                                * ( x [ 6 ]
	                                                                        - _params [ 7 ]
	                                                                                * _params [ 1 ]
	                                                                                * _params [ 2 ]
	                                                                                / pow (
	                                                                                        _params [ 0 ],
	                                                                                        2.0 ) ) )
	                                        + ( _params [ 10 ] - _params [ 11 ]
	                                                * _params [ 12 ] ) * v ) );

	dx [ 15 ] = _params [ 13 ] * x [ 9 ] + _params [ 14 ] * x [ 11 ]
	        - _params [ 13 ] * x [ 13 ] - ( _params [ 14 ] + _params [ 15 ] )
	        * x [ 15 ] + 1.0 / _params [ 4 ] * u [ 0 ];
}

double yane::ExampleModels2::ARP::continuousCostfunction ( double t,
        double * x, double * u )
{
	return pow ( x [ 8 ] - x [ 9 ], 2.0 );
}

void yane::ExampleModels2::ARP::restrictionFunction ( double t, double * x,
        double * u, double * fx )
{
}

double yane::ExampleModels2::ARP::discreteCostfunction ( int length,
        int horizont, double * t, double * x, double * u )
{
	return 0.0;
}

void yane::ExampleModels2::ARP::getObjectiveWeight (
        double & weight_continuousCostfunction,
        double & weight_discreteCostfunction )
{
	weight_continuousCostfunction = 1.0;
	weight_discreteCostfunction = 0.0;
}

void yane::ExampleModels2::ARP::getDefaultState ( double * x )
{
	x [ 0 ] = 0.0;
	x [ 1 ] = 0.0;
	x [ 2 ] = 0.0;
	x [ 3 ] = 0.0;
	x [ 4 ] = 0.0;
	x [ 5 ] = 0.0;
	x [ 6 ] = 0.0;
	x [ 7 ] = 0.0;
	x [ 8 ] = 10.0;
	x [ 9 ] = 10.0;
	x [ 10 ] = 0.0;
	x [ 11 ] = 0.0;
	x [ 12 ] = 0.0;
	x [ 13 ] = 0.0;
	x [ 14 ] = 0.0;
	x [ 15 ] = 0.0;
}

void yane::ExampleModels2::ARP::getDefaultControl ( double * u )
{
	u [ 0 ] = 0.0;
}

void yane::ExampleModels2::ARP::eventBeforeMPC ( int horizont, double * t,
        double * x, double * sdatavalues )
{
	// 	_params[22] = 0.0;					// Tracking zero solution
	// 	_params[23] = 0.0;
	// 	_params[24] = 0.0;
	// 	_params[25] = 0.0;
	// 	_params[26] = 0.0;

	if ( t [ 0 ] < 5.0 ) // Tracking reference
	{
		_params [ 22 ] = 10.0;
		_params [ 23 ] = 0.0;
		_params [ 24 ] = 0.0;
		_params [ 25 ] = 0.0;
		_params [ 26 ] = 0.0;
	}
	else if ( t [ 0 ] < 9.0 )
	{
		_params [ 22 ] = 0.0;
		_params [ 23 ] = 0.0;
		_params [ 24 ] = 0.0;
		_params [ 25 ] = 0.0;
		_params [ 26 ] = 0.0;
	}
	else if ( t [ 0 ] < 10.0 )
	{
		_params [ 22 ] = 10.0;
		_params [ 23 ] = 0.0;
		_params [ 24 ] = 0.0;
		_params [ 25 ] = 0.0;
		_params [ 26 ] = 0.0;
	}
	else
	{
		_params [ 22 ] = 0.0;
		_params [ 23 ] = 0.0;
		_params [ 24 ] = 0.0;
		_params [ 25 ] = 0.0;
		_params [ 26 ] = 0.0;
	}
}

void yane::ExampleModels2::ARP::setModelParameter ( void )
{
	setDoubleParameter ( 0, 5.0 ); // Total mass
	setDoubleParameter ( 1, 0.5 ); // Mass of arm
	setDoubleParameter ( 2, 0.3 ); // Radius
	setDoubleParameter ( 3, 0.06 ); // Moment of arm
	setDoubleParameter ( 4, 0.005 ); // Moment of rotor
	setDoubleParameter ( 5, 64.0 ); // Spring constant k1
	setDoubleParameter ( 6, 8.0 ); // Spring constant k3

	setDoubleParameter ( 7, 1.60 );
	setDoubleParameter ( 8, 0.04 );
	setDoubleParameter ( 9, 0.007 );

	setDoubleParameter ( 10, _params [ 6 ] * _params [ 0 ] / ( _params [ 0 ]
	        * _params [ 3 ] - pow ( _params [ 1 ] * _params [ 2 ], 2.0 ) ) );
	setDoubleParameter ( 11, ( _params [ 8 ] * pow ( _params [ 0 ], 2.0 )
	        + _params [ 7 ] * pow ( _params [ 1 ] * _params [ 2 ], 2.0 ) )
	        / ( _params [ 0 ] * ( _params [ 0 ] * _params [ 3 ] - pow (
	                _params [ 1 ] * _params [ 2 ], 2.0 ) ) ) );
	setDoubleParameter ( 12, _params [ 8 ] * _params [ 0 ] / ( _params [ 0 ]
	        * _params [ 3 ] - pow ( _params [ 1 ] * _params [ 2 ], 2.0 ) ) );
	setDoubleParameter ( 13, _params [ 6 ] / _params [ 4 ] );
	setDoubleParameter ( 14, _params [ 8 ] / _params [ 4 ] );
	setDoubleParameter ( 15, _params [ 9 ] / _params [ 4 ] );
	setDoubleParameter ( 16, _params [ 1 ] * _params [ 2 ] / ( _params [ 0 ]
	        * _params [ 3 ] - pow ( _params [ 1 ] * _params [ 2 ], 2.0 ) )
	        * _params [ 5 ] );
	setDoubleParameter ( 17, _params [ 1 ] * _params [ 2 ] / ( _params [ 0 ]
	        * _params [ 3 ] - pow ( _params [ 1 ] * _params [ 2 ], 2.0 ) )
	        * _params [ 7 ] );

	setDoubleParameter ( 18, 10000.0 ); // Design parameter 1 for continous control
	setDoubleParameter ( 19, 3500.0 ); // Design parameter 2 for continous control
	setDoubleParameter ( 20, 500.0 ); // Design parameter 3 for continous control
	setDoubleParameter ( 21, 35.0 ); // Design parameter 4 for continous control

	setDoubleParameter ( 22, 10.0 ); // Initial values for tracking
	setDoubleParameter ( 23, 0.0 );
	setDoubleParameter ( 24, 0.0 );
	setDoubleParameter ( 25, 0.0 );
	setDoubleParameter ( 26, 0.0 );
}

const yane::Utils::Uuid & yane::ExampleModels2::ARP::id ( )
{
	return ID;
}
