#ifndef _MAYBE_H__
#define _MAYBE_H__

#pragma once
#include <vector>
#include <iostream>
#include "../Circuit/InvalidOperationException.h"

namespace Utility
{
	template<typename T>
	class IMaybe
	{
	public:
		virtual T getValue() const = 0;
	};

	template<typename T>
	class Nothing : public IMaybe<T>
	{
	public:
		T getValue() const
		{
			throw new InvalidOperationException();
		}
	};

	template<typename T>
	class Just : public IMaybe<T>
	{
	public:
		const T val;

		Just(T val) : val(val)
		{
		}

		T getValue() const
		{
			return val;
		}
	};
}

#endif
