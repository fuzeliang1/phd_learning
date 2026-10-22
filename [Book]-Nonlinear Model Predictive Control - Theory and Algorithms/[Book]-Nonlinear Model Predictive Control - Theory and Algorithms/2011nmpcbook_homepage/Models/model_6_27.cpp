/***************************************************************************
 *                                                                         *
 * Advection-Diffusion PDE                                                 *
 *                                                                         *
 * Authors:                                                                *
 *   Juergen Pannek <juergen.pannek@uni-bayreuth.de>                       *
 *                                                                         *
 ***************************************************************************/

/**
 @file
 @author Juergen Pannek <juergen.pannek@uni-bayreuth.de>
 */

#include "model_6_27.h"

#include <yane/odesolve.h>
#include <iostream>
#include <cmath>

using namespace std;

#define PI 3.141592654
#define INF 1.0E19
#define LENGTH 1.0

yane::ExampleModels2::SemiHeat::SemiHeat ( int dimensionState,
        int dimensionControl ) :
	yane::Model::Model ( new yane::OdeSolve::DoPri853 ( ),
	        new yane::OdeSolve::DoPriConfig ( ), dimensionState,
	        dimensionControl, 0, 0 )
{
	odeConfig ( )->setTolerance ( 1E-10, 1E-10 );

	if ( dimensionControl == 2 )
	{
		_distributedControl = false;
	}
	else if ( dimensionControl == dimensionState )
	{
		_distributedControl = true;
	}
	else
	{
		cout
		        << "Error on construction object, control dimension must be either 2 or equal to state dimension"
		        << endl;
	}

	_h1norm = false;
}

yane::ExampleModels2::SemiHeat::~SemiHeat ( )
{
	delete odeSolver ( );
	delete odeConfig ( );
}

void yane::ExampleModels2::SemiHeat::odeFunction ( double t, double * x,
        double * u, double * dx )
{
	double theta = 0.1;
	double rho = 10.0;
	double grid = LENGTH / double ( _dimension + 1 );

	if ( ! _distributedControl )
	{
		dx [ 0 ] = theta * ( u [ 0 ] - 2.0 * x [ 0 ] + x [ 1 ] ) / pow ( grid,
		        2 ) - ( x [ 1 ] - x [ 0 ] ) / grid + rho * semiFunction (
		        x [ 0 ] );
		for ( int i = 1; i < _dimension - 1; i++ )
		{
			dx [ i ] = theta * ( x [ i - 1 ] - 2.0 * x [ i ] + x [ i + 1 ] )
			        / pow ( grid, 2 ) - ( x [ i + 1 ] - x [ i ] ) / grid + rho
			        * semiFunction ( x [ i ] );
		}
		dx [ _dimension - 1 ] = theta * ( x [ _dimension - 2 ] - 2.0
		        * x [ _dimension - 1 ] + u [ 1 ] ) / pow ( grid, 2 ) - ( u [ 1 ]
		        - x [ _dimension - 1 ] ) + rho * semiFunction ( x [ _dimension
		        - 1 ] );
	}
	else
	{
		dx [ 0 ] = theta * ( - 2.0 * x [ 0 ] + x [ 1 ] ) / pow ( grid, 2 )
		        + ( x [ 1 ] - x [ 0 ] ) / grid + rho * semiFunction ( x [ 0 ] )
		        + u [ 0 ];
		for ( int i = 1; i < _dimension - 1; i++ )
		{
			dx [ i ] = theta * ( x [ i - 1 ] - 2.0 * x [ i ] + x [ i + 1 ] )
			        / pow ( grid, 2 ) + ( x [ i + 1 ] - x [ i ] ) / grid + rho
			        * semiFunction ( x [ i ] ) + u [ i ];
		}
		dx [ _dimension - 1 ] = theta * ( x [ _dimension - 2 ] - 2.0
		        * x [ _dimension - 1 ] ) / pow ( grid, 2 ) + ( - x [ _dimension
		        - 1 ] ) + rho * semiFunction ( x [ _dimension - 1 ] )
		        + u [ _dimension - 1 ];
	}
}

double yane::ExampleModels2::SemiHeat::continuousCostfunction ( double t,
        double * x, double * u )
{
	double omega = 0.0;
	if ( _h1norm )
	{
		omega = 1.0;
	}
	double lambda = 0.1;
	double x_sum = 0.0;
	double dx_sum = 0.;
	double u_sum = 0.0;
	double dx = LENGTH / double ( _dimension + 1 );
	double fx = 0.0;

	if ( ! _distributedControl )
	{
		for ( int j = 0; j < ( _dimension ); j++ )
		{
			x_sum += pow ( x [ j ], 2 );
		}
		x_sum += 0.5 * ( pow ( u [ 0 ], 2 ) + pow ( u [ 1 ], 2 ) );

		for ( int j = 0; j < ( _dimension - 1 ); j++ )
		{
			dx_sum += pow ( ( x [ j + 1 ] - x [ j ] ) / dx, 2 );
		}
		dx_sum += 0.5 * ( pow ( ( x [ 0 ] - u [ 0 ] ) / dx, 2 ) + pow (
		        ( u [ 1 ] - x [ _dimension - 1 ] ) / dx, 2 ) );

		for ( int j = 0; j < ( _ctrl_dimension ); j++ )
		{
			u_sum += 0.5 * pow ( u [ j ], 2 );
		}
		fx = ( x_sum + omega * dx_sum ) * dx + lambda * u_sum;
		return fx;
	}
	else
	{
		for ( int j = 0; j < ( _dimension ); j++ )
		{
			x_sum += pow ( x [ j ], 2 );
		}

		for ( int j = 0; j < ( _dimension - 1 ); j++ )
		{
			dx_sum += pow ( ( x [ j + 1 ] - x [ j ] ) / dx, 2 );
		}
		dx_sum += pow ( x [ 0 ] / dx, 2 ) + pow ( x [ _dimension - 1 ] / dx, 2 );

		for ( int j = 0; j < ( _ctrl_dimension ); j++ )
		{
			u_sum += pow ( u [ j ], 2 );
		}
		fx = ( x_sum + omega * dx_sum + lambda * u_sum ) * dx;
		return fx;
	}
}

void yane::ExampleModels2::SemiHeat::restrictionFunction ( double t,
        double * x, double * u, double * fx )
{
}

double yane::ExampleModels2::SemiHeat::discreteCostfunction ( int length,
        int horizon, double * t, double * x, double * u )
{
	double omega = 0.0;
	if ( _h1norm )
	{
		omega = 1.0;
	}
	double lambda = 0.1;
	double x_sum = 0.0;
	double x_norm = 0.0;
	double dx_sum = 0.0;
	double dx_norm = 0.0;
	double u_sum = 0.0;
	double u_norm = 0.0;
	double dx = LENGTH / double ( _dimension + 1 );
	double fx = 0.0;

	if ( ! _distributedControl )
	{
		for ( int i = 0; i < horizon; i++ )
		{
			for ( int j = 0; j < ( _dimension ); j++ )
			{
				x_sum += pow ( x [ i * _dimension + j ], 2 );
			}
			x_sum += 0.5 * ( pow ( u [ i * _ctrl_dimension + 0 ], 2 ) + pow (
			        u [ i * _ctrl_dimension + 1 ], 2 ) );

			for ( int j = 0; j < ( _dimension - 1 ); j++ )
			{
				dx_sum += pow ( ( x [ i * _dimension + j + 1 ] - x [ i
				        * _dimension + j ] ) / dx, 2 );
			}
			dx_sum += 0.5 * ( pow (
			        ( x [ i * _dimension + 0 ] - u [ 0 ] ) / dx, 2 )
			        + pow ( ( u [ 1 ] - x [ i * _dimension + _dimension - 1 ] )
			                / dx, 2 ) );

			for ( int j = 0; j < ( _ctrl_dimension ); j++ )
			{
				u_sum += 0.5 * pow ( u [ i * _ctrl_dimension + j ], 2 );
			}
		}
		fx = ( x_sum + omega * dx_sum ) * dx + lambda * u_sum;
		return fx;
	}
	else
	{
		for ( int i = 0; i < horizon; i++ )
		{
			for ( int j = 0; j < ( _dimension ); j++ )
			{
				x_sum += pow ( x [ i * _dimension + j ], 2 );
			}

			for ( int j = 0; j < ( _dimension - 1 ); j++ )
			{
				dx_sum += pow ( ( x [ i * _dimension + j + 1 ] - x [ i
				        * _dimension + j ] ) / dx, 2 );
			}
			dx_sum += pow ( x [ i * _dimension + 0 ] / dx, 2 ) + pow ( x [ i
			        * _dimension + _dimension - 1 ] / dx, 2 );

			for ( int j = 0; j < ( _ctrl_dimension ); j++ )
			{
				u_sum += pow ( u [ i * _ctrl_dimension + j ], 2 );
			}
		}

		fx = ( x_sum + omega * dx_sum + lambda * u_sum ) * dx;
		return fx;
	}
}

void yane::ExampleModels2::SemiHeat::getObjectiveWeight (
        double & weight_continuousCostfunction,
        double & weight_discreteCostfunction )
{
	weight_continuousCostfunction = 1.0;
	weight_discreteCostfunction = 0.0;
}

void yane::ExampleModels2::SemiHeat::getDefaultState ( double * x )
{
	int i;
	for ( i = 0; i < _dimension; i++ )
	{
		x [ i ] = 0.1;
	}
	for ( i = 0; i < floor ( double ( _dimension ) / 3.3 ); i++ )
	{
		x [ i ] = - 0.1;
	}
}

void yane::ExampleModels2::SemiHeat::getDefaultControl ( double * u )
{
	for ( int i = 0; i < _ctrl_dimension; i++ )
	{
		u [ i ] = 0;
	}
}

double yane::ExampleModels2::SemiHeat::semiFunction ( double x )
{
	return ( x - pow ( x, 3 ) );
}

void yane::ExampleModels2::SemiHeat::eventBeforeMPC ( int horizon,
        double * timesteps, double * x, double * sdatavalues )
{
}

