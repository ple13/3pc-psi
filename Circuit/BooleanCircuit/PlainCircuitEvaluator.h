#ifndef PLAIN_CIRCUIT_EVALUATOR_H__
#define PLAIN_CIRCUIT_EVALUATOR_H__

#pragma once

#include "InvalidOperationException.h"
#include "BooleanCircuit.h"
#include <boost/dynamic_bitset.hpp>

namespace Circuit
{
	class PlainBooleanCircuitEvaluator
	{
	public:
		boost::dynamic_bitset<> *EvaluateBooleanCircuit(BooleanCircuit *circuit, boost::dynamic_bitset<> *input);
	};
}

#endif
