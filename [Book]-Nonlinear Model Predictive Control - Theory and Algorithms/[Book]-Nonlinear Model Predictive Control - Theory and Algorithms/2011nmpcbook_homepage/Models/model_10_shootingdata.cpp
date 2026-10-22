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

#include "model_10_shootingdata.h"
#include "model_10.h"
#include <iostream>
#include <yane/model.h>
#include <yane/utils.h>

yane::ExampleModels2::InvertedPendulumShootingData::InvertedPendulumShootingData (
        yane::Model::Model * pendulum )
{
	_pendulum = pendulum;
	_numberOfShootingNodes = 0;
}

void yane::ExampleModels2::InvertedPendulumShootingData::applyShootingData (
        int horizon, double * t, double * x, double * sdatavalues )
{
	for ( int i = 0; i < _numberOfShootingNodes; i++ )
	{
		sdatavalues [ i ] = _shootingValue [ i ];
	}
}

int yane::ExampleModels2::InvertedPendulumShootingData::shootingDataLength (
        int horizon )
{
	return _numberOfShootingNodes;
}

int yane::ExampleModels2::InvertedPendulumShootingData::maxShootingDataLength (
        int maxhorizon )
{
	return _numberOfShootingNodes;
}

void yane::ExampleModels2::InvertedPendulumShootingData::shootingDataInfo (
        int horizon, yane::MPC::T_SHOOTINGDATA * sdata )
{
	for ( int i = 0; i < _numberOfShootingNodes; i++ )
	{
		sdata [ i ].horizonindex = _shootingHorizonIndex [ i ];
		sdata [ i ].varindex = _shootingVarIndex [ i ];
	}
}

yane::ExampleModels2::InvertedPendulumShootingData::~InvertedPendulumShootingData ( )
{
}

void yane::ExampleModels2::InvertedPendulumShootingData::setShooting (
        int horizon, int numberOfShootingNodes, int * shootingHorizonIndex,
        int * shootingVarIndex, double * shootingValue )
{
	_numberOfShootingNodes = numberOfShootingNodes;
	for ( int i = 0; i < _numberOfShootingNodes; i++ )
	{
		if ( shootingHorizonIndex [ i ] > horizon - 1
		        || shootingHorizonIndex [ i ] < 1 )
		{
			throw yane::Utils::ValueException (
			        __DEBUGMESS (
			                "Index of shooting node on horizon must be between 1 and the length of the horizon" ) );
		}
		if ( shootingVarIndex [ i ] >= _pendulum->dimensionState ( )
		        || shootingVarIndex < 0 )
		{
			throw yane::Utils::ValueException (
			        __DEBUGMESS (
			                "Index of shooting node on state dimension must be between 0 and the dimension of the state" ) );
		}
	}
	_shootingHorizonIndex = shootingHorizonIndex;
	_shootingVarIndex = shootingVarIndex;
	_shootingValue = shootingValue;
}

