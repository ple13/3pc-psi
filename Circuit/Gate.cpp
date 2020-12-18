#include "Gate.h"
#include <assert.h>
#include <vector>
#include <string>

namespace Circuit
{

	UnitaryGate::UnitaryGate(int inputWire, int outputWire) : InputWire(inputWire), OutputWire(outputWire){}
	
	UnitaryGate::UnitaryGate(const UnitaryGate &obj)
	{
		this->InputWire = obj.InputWire;
		this->OutputWire = obj.OutputWire;
	}
	
	UnitaryGate::UnitaryGate(const IGate &other)
	{
		const UnitaryGate &obj = dynamic_cast<const UnitaryGate &>(other);
		*this = obj;
	}
	
	std::string UnitaryGate::ToString()
	{
		return std::to_string(InputWire) + "|" + std::to_string(OutputWire);
	}

	BinaryGate::BinaryGate(int LeftWire, int RightWire, int OutputWire) : LeftWire(LeftWire), RightWire(RightWire), OutputWire(OutputWire){}
	
	BinaryGate::BinaryGate(const BinaryGate &obj)
	{
		this->LeftWire = obj.LeftWire;
		this->RightWire = obj.RightWire;
		this->OutputWire = obj.OutputWire;
	}
	
	BinaryGate::BinaryGate(const IGate &other)
	{
		const BinaryGate &obj = dynamic_cast<const BinaryGate &>(other);
		this->LeftWire = obj.LeftWire;
		this->RightWire = obj.RightWire;
		this->OutputWire = obj.OutputWire;
	}
	
	std::string BinaryGate::ToString()
	{
		return std::to_string(LeftWire) + "|" + std::to_string(RightWire) + "|" + std::to_string(OutputWire);
	}
	
	NaryGate::NaryGate(std::vector<int> inputWires, int outputWire) : InputWires(std::move(inputWires)), OutputWire(outputWire){}
	
	std::string NaryGate::ToString()
	{
		std::string str;
		for(int idx = 0; idx < InputWires.size(); idx++)
		{
			str = str + std::to_string(InputWires[idx]) + "|";
		}
		str = str + std::to_string(OutputWire);
		return str;
	}
}
