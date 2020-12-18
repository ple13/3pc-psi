#include "BooleanGate.h"


namespace Circuit
{

	AndGate::AndGate(int LeftWire, int RightWire, int OutputWire) : BinaryGate(LeftWire, RightWire, OutputWire)
	{
	}

	XorGate::XorGate(int LeftWire, int RightWire, int OutputWire) : BinaryGate(LeftWire, RightWire, OutputWire)
	{
	}

	NotGate::NotGate(int inputWire, int outputWire) : UnitaryGate(inputWire, outputWire)
	{
	}
}
