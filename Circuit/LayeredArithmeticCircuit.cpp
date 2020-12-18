#include <algorithm>
#include "LayeredArithmeticCircuit.h"

namespace Circuit
{

	ArithmeticLayer::ArithmeticLayer(std::vector<AdditionGate *> AddGates, 
					 std::vector<SubtractionGate *> SubGates,
				         std::vector<MultiplicationGate *> MulGates, 
					 std::vector<DotProductGate *> DotGates,
				         std::vector<DivisionGate *> DivGates,
					 std::vector<NaryAdditionGate *> NAddGates,
					 std::vector<NaryDotGate *> NDotGates,
				         std::vector<ConstantAdditionGate *> CAddGates, 
				         std::vector<ConstantMultiplicationGate *> CMulGates) : 
				         AddGates(std::move(AddGates)), SubGates(std::move(SubGates)), 
				         MulGates(std::move(MulGates)), DotGates(std::move(DotGates)),
				         DivGates(std::move(DivGates)), NDotGates(std::move(NDotGates)), 
				         NAddGates(std::move(NAddGates)), CAddGates(std::move(CAddGates)),
				         CMulGates(std::move(CMulGates))
	{
	}

	std::string ArithmeticLayer::ToString()
	{
		std::string builder;
		
		for (auto &add : AddGates)
		{
			builder = builder + "ADD" + "|" + add->ToString() + "\n";
		}
		
		for (auto &sub : SubGates)
		{
			builder = builder + "SUB" + "|" + sub->ToString() + "\n";
		}

		for  (auto &mul : MulGates)
		{
			builder = builder + "MUL" + "|" + mul->ToString() + "\n";
		}
		
		for  (auto &dot : DotGates)
		{
			builder = builder + "DOT" + "|" + dot->ToString() + "\n";
		}
		
		for (auto &div : DivGates)
		{
			builder = builder + "DIV" + "|" + div->ToString() + "\n";
		}
		
		for (auto &nadd : NAddGates)
		{
			builder = builder + "NADD" + "|" + nadd->ToString() + "\n";
		}
		
		for (auto &ndot : NDotGates)
		{
			builder = builder + "NDOT" + "|" + ndot->ToString() + "\n";
		}

		for  (auto &cadd : CAddGates)
		{
			builder = builder + "CADD" + "|" + cadd->ToString() + "\n";
		}

		for  (auto &cmul : CMulGates)
		{
			builder = builder + "CMUL" + "|" + cmul->ToString() + "\n";
		}
		
		return builder;
	}

	LayeredArithmeticCircuit::LayeredArithmeticCircuit(std::vector<ArithmeticLayer *> &layers, int depth, int numWires, int numMulGates, int inputLength, int outputLength) : LayeredCircuit<ArithmeticLayer>(layers, depth, numWires, numMulGates, inputLength, outputLength)
	{
	}

	std::string LayeredArithmeticCircuit::ToString()
	{
		std::string builder;
		
		builder += "\n";//		f("");
		builder += "---------------------------------------\n"; //f("---------------------------------------");

		for (int i = 0; i < Depth; i++)
		{
			ArithmeticLayer *layer = this->operator[](i);
			builder = builder + "Layer " + std::to_string(i) + "\n"; //f("Layer " + std::to_wstring(i));
			builder = builder + layer->ToString() + "\n"; //builder->append(layer->ToString());
			builder += "---------------------------------------\n"; //f("---------------------------------------");
		}

		return builder;
	}
}
