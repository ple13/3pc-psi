#ifndef STRING_HELPER_H__
#define STRING_HELPER_H__

#pragma once

#include <string>
#include <sstream>
#include <vector>
#include <iostream>

class StringHelper
{
public:
	static std::vector<std::string> split(const std::string &source, char delimiter)
	{
		std::vector<std::string> output;
		std::istringstream ss(source);
		std::string nextItem;

		while (std::getline(ss, nextItem, delimiter))
		{
			output.push_back(nextItem);
		}

		return output;
	}
};

#endif
