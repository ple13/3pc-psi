#ifndef STRING_HELPER_H__
#define STRING_HELPER_H__

#pragma once

#include <string>
#include <sstream>
#include <vector>
#include <iostream>

std::ostream&
operator<<( std::ostream& dest, __uint128_t value )
{
    std::ostream::sentry s( dest );
    if ( s ) {
        __uint128_t tmp = value < 0 ? -value : value;
        char buffer[ 128 ];
        char* d = std::end( buffer );
        do
        {
            -- d;
            *d = "0123456789"[ tmp % 10 ];
            tmp /= 10;
        } while ( tmp != 0 );
        if ( value < 0 ) {
            -- d;
            *d = '-';
        }
        int len = std::end( buffer ) - d;
        if ( dest.rdbuf()->sputn( d, len ) != len ) {
            dest.setstate( std::ios_base::badbit );
        }
    }
    return dest;
}

void compactData(unsigned char *out, const std::vector<__uint128_t>& in, int nBytes)
{
	// Compact the bytes
	unsigned char *ptr = 0;
	int count = 0;
	
	ptr = (unsigned char *)(in.data());
	
	for(int idx = 0; idx < in.size(); idx++)
	{
		for(int kdx = 0; kdx < nBytes; kdx++)
		{
			out[count] = *ptr;
			count++;
			ptr++;
		}
		ptr += (16 - nBytes);
	}
}

void uncompactData(std::vector<__uint128_t>& out, unsigned char *in, int nBytes)
{
	// Compact the bytes
	unsigned char *ptr = 0;
	int count = 0;
	
	ptr = (unsigned char *)(out.data());
	
	for(int idx = 0; idx < out.size(); idx++)
	{
		for(int kdx = 0; kdx < nBytes; kdx++)
		{
			*ptr = in[count];
			count++;
			ptr++;
		}
		for(int kdx = nBytes; kdx < 16; kdx++)
		{
			*ptr = 0; 
			ptr++;
		}
	}
}

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
