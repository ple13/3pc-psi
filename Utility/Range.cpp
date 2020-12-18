#include <cassert>
#include "Range.h"

namespace Utility
{

	Range::Range(int start, int length) : Start(start), Length(length)
	{
		assert(start >= 0);
		assert(length > 0);
	}
}
