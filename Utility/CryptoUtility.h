#ifndef CRYPTO_UTILITY_H__
#define CRYPTO_UTILITY_H__

#pragma once
#include <cstddef>
#include <cassert>
#include <stdint.h>
#include <vector>
#include <bitset>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

#include "Communicator.h"
#include "../GraphData.hpp"

using namespace NTL;


#define PRECISION_BIT_LENGTH 16

namespace Utility
{
	class VectorOperation
	{
	public:
		static std::vector<int64_t> Add(const std::vector<int64_t>& x, const std::vector<int64_t>& y);
		static std::vector<int64_t> Sub(const std::vector<int64_t>& x, const std::vector<int64_t>& y);
		static std::vector<int64_t> Mul(const int64_t x, const std::vector<int64_t>& y);
		static std::vector<int64_t> Mul(const std::vector<int64_t>& x, const int64_t y);
		static std::vector<int64_t> Mul(const std::vector<int64_t>& x, const std::vector<int64_t>& y);
		static std::vector<int64_t> Dot(const std::vector<int64_t>& x, const std::vector<int64_t>& y);
	};
	
	class CryptoUtility
	{
	public:
		static std::vector<int64_t> SampleSmallInput(int length);
		
		static std::vector<unsigned char> SampleByteArray(int length);

		static std::vector<unsigned char> ComputeHash(std::vector<unsigned char> msg);
		
		static std::vector<unsigned char> ComputeHash(unsigned char *string, int size);
		
		static int32_t * SampleInt32Array(int length);
		
		static std::vector<int64_t>  SampleInt64Array(int length);
		
		static std::vector<int> buildMaskIndex(const std::vector<int>& itemsPerUser);
		
		template<typename T>
		static void shuffle(std::vector<T>& data, const std::vector<unsigned char>& seed);

		template<typename T>
		static void shuffle(std::vector<T>& data, __int128& seed);
	};
	
	class ArithmeticOperation
	{
	public:
		static int64_t add(int64_t x, int64_t y){return x + y;}
		static int64_t sub(int64_t x, int64_t y){return x - y;}
		static int64_t mul(int64_t x, int64_t y){return (x * y) >> PRECISION_BIT_LENGTH;}
		static int64_t div(int64_t x, int64_t y){return (x << PRECISION_BIT_LENGTH) / y;}
	};
	
	class ArrayEncoder
	{
	public:
		static std::vector<int64_t> UCharVec2Int64tVec(std::vector<unsigned char> array);
		
		static std::vector<unsigned char> Encode(std::vector<int64_t> array);		
		static std::vector<unsigned char> Encode(std::vector<uint64_t> array);
		
		static std::vector<unsigned char> Hash(std::vector<int64_t> array);
		static std::vector<unsigned char> Hash(std::vector<uint64_t> array);
		
		static std::vector<int64_t>  Decodeint64_t(std::vector<unsigned char> &array);
		
		static std::vector<unsigned char> EncodeInt64Array(std::vector<std::vector<int64_t> > array);
		static std::vector<unsigned char> EncodeUInt64Array(std::vector<std::vector<uint64_t> > array);
		
		static std::vector<std::vector<int64_t> > DecodeInt64Array(std::vector<unsigned char> array);
		static std::vector<std::vector<uint64_t> > DecodeUInt64Array(std::vector<unsigned char> array);
		
		/// Convert a 4-byte array to unsigned int value
		static int byteArrayToInt32(unsigned char *byteArray);

		/// Convert an unsigned int to 4-byte array
		static unsigned char * int32ToByteArray(unsigned int val);
		
		static std::vector<unsigned char> MacVectorToByteArray(std::vector<uint64_t> MACedEdges);
		
		static std::vector<unsigned char> MacVectorToByteArray(uint64_t *MACedEdges, int size);
		
		static std::vector<uint64_t> ByteArrayToMacVector(std::vector<unsigned char> byteArray);
		
		static void ZZArray2ByteArray(unsigned char *dest, const ZZ *src,    int nBytes, int size);
		static void ZZArray2ByteArray(unsigned char *dest, const ZZ_p *src,  int nBytes, int size);
		
		static void ZZpX2ByteArray(unsigned char *dest, const ZZ_pX& src, int nBytes, int size);
		static void ByteArray2ZZpX(ZZ_pX& dest, unsigned char *src, int nBytes, int size);
		
		static void ByteArray2ZZArray(ZZ *dest,   unsigned char *src, int nBytes, int size);
		static void ByteArray2ZZArray(ZZ_p *dest, unsigned char *src, int nBytes, int size);
	};
	
	// 40:  	87, 167, 195, 203, 213, 285, 293, 299, 389, 437
	class DistributedMACGenerator
	{
	public:
		DistributedMACGenerator(uint64_t alpha, uint64_t beta);
		
		uint64_t computeMAC(int party, std::vector<uint64_t> &data, uint64_t random);
		
		__int128 Z64;
	private:
		const uint64_t p = 0x10000000000 - 293;
		
		uint64_t alpha;
		uint64_t beta;
		
		std::vector<uint64_t> coeff;
	};
	
	class SPDZ2kMAC
	{
	public:
		SPDZ2kMAC(uint64_t alpha, Communicator *com);
		uint64_t computeMAC(int party, uint64_t data, uint64_t random);
		uint64_t singleCheck(int party, uint64_t mac, uint64_t share);
// 		uint64_t multipleCheck(int party, std::vector<uint64_t> mac, std::vector<uint64_t> share);
// 	private:
		uint64_t alpha;
		Communicator *com;
	};
	
	class TestUtility
	{
	public:
		static void PrintByteArray(const std::vector<unsigned char>& array, const std::string& str);
		static void PrintVector(const std::vector<int16_t>& input, const std::string& str);
		static void PrintVector(const std::vector<uint16_t>& input, const std::string& str);
		static void PrintVector(const std::vector<int64_t>& input, const std::string& str);
		static void PrintVector(const std::vector<uint64_t>& input, const std::string& str);
		static void PrintVector(const std::vector<std::vector<int64_t> >& input, const std::string& str, float scale);
		static void PrintVector(const std::vector<std::vector<uint64_t> >& input, const std::string& str, float scale);
		static void PrintMask(const std::vector<int64_t>& masks, const std::vector<int>& maskIndex, const std::vector<int>& itemsPerUser, const std::string& str);
		static std::vector<std::vector<int64_t> > GenerateInput(const std::vector<int>& itemsPerUser, int inputLength);
	};
}

#endif
