#ifndef RANGE_H__
#define RANGE_H__

#pragma once

namespace Utility
{
	class IRange
	{
	public:
		virtual ~IRange(){};
	};

	class EmptyRange : public IRange
	{
	public:
		~EmptyRange(){}
	};

	class Range : public IRange
	{
	public:
		const int Start;
		const int Length;

		Range(int start, int length);
		~Range(){}
	};
}

#endif
