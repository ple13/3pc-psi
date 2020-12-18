#ifndef CIRCUIT_COMPUTATION_H__
#define CIRCUIT_COMPUTATION_H__

#include "Shares.h"
#include "Triplets.h"

Shares RippleCarryAdder(const Shares& x, const Shares& y)
{
	Shares s = x;
	Shares c = y;
	Shares temp1, temp2;
	
	for(int idx = 0; idx < 80; idx++)
	{
		temp1 = Shares::ComputeXOR(c, s);
		temp2 = Shares::ComputeAND(c, s);
		
		s = temp1;
		c = temp2;
		
		c.LeftShiftOne();
	}
	
	return s;
}

#endif
