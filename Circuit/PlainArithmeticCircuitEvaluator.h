#ifndef PLAIN_ARITHMETIC_CIRCUIT_EVALUATOR_H__
#define PLAIN_ARITHMETIC_CIRCUIT_EVALUATOR_H__

#pragma once
#include <vector>
#include "ArithmeticCircuit.h"
#include "LayeredArithmeticCircuit.h"

namespace Circuit
{
	class PlainArithmeticCircuitEvaluator
	{
	public:
		void EvaluateArithmeticCircuit(ArithmeticCircuit *circuit, LayeredArithmeticCircuit *lc, const std::vector<int64_t>& input, std::vector<int64_t>& evaluation);
		
		void EvaluateArithmeticCircuit(ArithmeticCircuit *circuit, LayeredArithmeticCircuit *lc, const std::vector<std::vector<int64_t> >& input, const std::vector<int>& itemsPerUser, std::vector<std::vector<int64_t> >& evaluation);
	};
}

#endif
