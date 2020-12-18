#include <cstring>
#include <iostream>
#include <sstream>
#include <random>
#include <chrono>
#include <algorithm>
#include <iterator>
#include <openssl/rand.h>
#include <openssl/sha.h>
#include "CryptoUtility.h"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
using namespace NTL;

namespace Utility
{
	std::vector<int64_t> VectorOperation::Add(const std::vector<int64_t>& x, const std::vector<int64_t>& y)
	{
		std::vector<int64_t> ret(x.size());
		
		for(int idx = 0; idx < x.size(); idx++)
		{
			ret[idx] = x[idx] + y[idx];
		}
		
		return ret;
	}
	
	std::vector<int64_t> VectorOperation::Sub(const std::vector<int64_t>& x, const std::vector<int64_t>& y)
	{
		std::vector<int64_t> ret(x.size());
		
		for(int idx = 0; idx < x.size(); idx++)
		{
			ret[idx] = x[idx] - y[idx];
		}
		
		return ret;
	}
	
	std::vector<int64_t> VectorOperation::Mul(const int64_t x, const std::vector<int64_t>& y)
	{
		std::vector<int64_t> ret(y.size());
		
		for(int idx = 0; idx < y.size(); idx++)
		{
			ret[idx] = ArithmeticOperation::mul(x, y[idx]);
		}
		
		return ret;
	}
	
	std::vector<int64_t> VectorOperation::Mul(const std::vector<int64_t>& x, const int64_t y)
	{
		std::vector<int64_t> ret(x.size());
		
		for(int idx = 0; idx < x.size(); idx++)
		{
			ret[idx] = ArithmeticOperation::mul(y, x[idx]);
		}
		
		return ret;
	}
	
	std::vector<int64_t> VectorOperation::Mul(const std::vector<int64_t>& x, const std::vector<int64_t>& y)
	{
		std::vector<int64_t> ret(x.size());
		
		for(int idx = 0; idx < x.size(); idx++)
		{
			ret[idx] = ArithmeticOperation::mul(x[idx], y[idx]);
		}
		
		return ret;
	}
	
	std::vector<int64_t> VectorOperation::Dot(const std::vector<int64_t>& x, const std::vector<int64_t>& y)
	{
		std::vector<int64_t> ret(1);
		
		ret[0] = 0;
		
		for(int idx = 0; idx < x.size(); idx++)
		{
			ret[0] += ArithmeticOperation::mul(x[idx], y[idx]);
		}
		
		return ret;
	}
	
	std::vector<int64_t> CryptoUtility::SampleSmallInput(int length){
		unsigned seed =  std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
		std::uniform_int_distribution<int> distribution(-3, 3);

		std::vector<int64_t> input;
		for(int idx = 0; idx < length; idx++)
		{
			int64_t rng = ((distribution(generator)) << 17);
			input.push_back(rng);
		}
		
		return input;
	}
	
	std::vector<unsigned char> CryptoUtility::SampleByteArray(int length)
	{
		unsigned char *sample = new unsigned char[length];
		RAND_bytes(sample, length);
		std::vector<unsigned char> ret(sample, sample + length);
		return ret;
	}

	std::vector<unsigned char> CryptoUtility::ComputeHash(std::vector<unsigned char> string)
	{
		unsigned char *hash = new unsigned char[SHA256_DIGEST_LENGTH];
		SHA256_CTX sha256;
		SHA256_Init(&sha256);
		SHA256_Update(&sha256, string.data(), string.size());
		SHA256_Final(hash, &sha256);
		std::vector<unsigned char> ret(hash, hash + SHA256_DIGEST_LENGTH);
		return ret;
	}
	
	std::vector<unsigned char> CryptoUtility::ComputeHash(unsigned char *string, int size)
	{
		unsigned char *hash = new unsigned char[SHA256_DIGEST_LENGTH];
		SHA256_CTX sha256;
		SHA256_Init(&sha256);
		SHA256_Update(&sha256, string, size);
		SHA256_Final(hash, &sha256);
		std::vector<unsigned char> ret(hash, hash + SHA256_DIGEST_LENGTH);
		return ret;
	}
	
	int32_t * CryptoUtility::SampleInt32Array(int length)
	{
		unsigned char *sample = new unsigned char[length*sizeof(int32_t)];
		RAND_bytes(sample, length*sizeof(int32_t));
		return (int32_t *)(sample);
	}
	
	std::vector<int64_t> CryptoUtility::SampleInt64Array(int length)
	{
		// Sample values in the range of 32 bits
		unsigned char *sample = new unsigned char[length*sizeof(int32_t)];
		RAND_bytes(sample, length*sizeof(int32_t));
		int32_t *ptr = (int32_t *)sample;
		
		std::vector<int64_t> ret(length);
		
		for(int idx = 0; idx < length; idx++)
		{
			ret[idx] = ptr[idx];
		}
		
		delete [] sample;
		
		return ret;
	}
	
	std::vector<int> CryptoUtility::buildMaskIndex(const std::vector<int>& itemsPerUser){
		// Count number of required random values: some wires need 10 masks, some 1
		int numUsers = itemsPerUser.size();
		int numEdges = 0;
		for(int idx = 0; idx < itemsPerUser.size(); idx++)
		{
			numEdges += itemsPerUser[idx];
		}
		
		int numMasks = 6*dimension*numUsers + (2*dimension + 4)*numEdges;
		int numWires = 6*numUsers + 6*numEdges;
		
		std::vector<int> maskIndex(numMasks);
		
		int shift = 0;
		int count = 0; // increment for output wire of each gate
		
// 		std::cout << "Masks for Input" << std::endl;
		// Go through each user
		for(int idx = 0; idx < numUsers; idx++)
		{
			maskIndex[shift] = count; count += dimension; shift++;		// u_idx
			for(int kdx = 0; kdx < itemsPerUser[idx]; kdx++)
			{
				maskIndex[shift] = count; count += dimension; shift++;	// v_idx_kdx
				maskIndex[shift] = count; count++; shift++;		// r_idx_kdx
				maskIndex[shift] = count; count++; shift++;		// isReal_idx_kdx
			}
		}
		
		// shift = numUsers + 3*numEdges
		// count = dimension*numUsers + (dimension + 2)numEdges
		
// 		std::cout << "Masks for DOT" << std::endl;
		// Assign masks for output wire of DOT PRODUCT GATE
		// Go through each user
		for(int idx = 0; idx < numUsers; idx++)
		{
			for(int kdx = 0; kdx < itemsPerUser[idx]; kdx++)
			{
				maskIndex[shift] = count; count++; shift++;		// DOT_idx_kdx
			}
		}
		
		// shift = numUsers + 4*numEdges
		// count = dimension*numUsers + (dimension + 3)numEdges
		
// 		std::cout << "Masks for SUB" << std::endl;
		// Assign masks for output wire of SUB GATE
		// Go through each user
		for(int idx = 0; idx < numUsers; idx++)
		{
			for(int kdx = 0; kdx < itemsPerUser[idx]; kdx++)
			{
				maskIndex[shift] = count; count++; shift++;		// SUB_idx_kdx
			}
		}
		
		// shift = numUsers + 5*numEdges
		// count = dimension*numUsers + (dimension + 4)numEdges
		
// 		std::cout << "Masks for MUL" << std::endl;
		// Assign masks for output wire of MUL GATE
		// Go through each user
		for(int idx = 0; idx < numUsers; idx++)
		{
			for(int kdx = 0; kdx < itemsPerUser[idx]; kdx++)
			{
				maskIndex[shift] = count; count += dimension; shift++;	// MUL_idx_kdx
			}
		}
		
		// shift = numUsers + 6*numEdges
		// count = dimension*numUsers + (2*dimension + 4)numEdges
		
// 		std::cout << "Masks for NDOT" << std::endl;
		// Assign masks for output wire of MUL GATE
		// Go through each user
		for(int idx = 0; idx < numUsers; idx++)
		{
			maskIndex[shift] = count; count += dimension; shift++;	// NDOT_idx
		}
		
		// shift = 2*numUsers + 5*numEdges
		// count = 2*dimension*numUsers + (2*dimension + 4)numEdges
		
// 		std::cout << "Masks for ADD" << std::endl;
		// Assign masks for output wire of NADD GATE, 2 CMUL GATEs, and 2 ADD GATEs
		// Go through each user

		for(int idx = 0; idx < 5*numUsers; idx++)
		{
			maskIndex[shift] = count; count += dimension; shift++;		// ADD_idx_kdx
		}
		
// 		std::cout << "count: " << count << std::endl;
		// shift = 6*numUsers + 6*numEdges
		// count = 7*dimension*numUsers + (2*dimension + 4)numEdges
		
		return maskIndex;
	}
	
	template<typename T>
	void CryptoUtility::shuffle(std::vector<T>& data, const std::vector<unsigned char>& seed)
	{
		__int128 *val = (__int128 *)seed.data();
		std::mt19937 g(*val);
		std::shuffle(data.begin(), data.end(), g);
	}

	template<typename T>
	void CryptoUtility::shuffle(std::vector<T>& data, __int128 & seed)
	{
		std::mt19937 g(seed);
		std::shuffle(data.begin(), data.end(), g);
	}
	
	std::vector<int64_t> ArrayEncoder::UCharVec2Int64tVec(std::vector<unsigned char> array)
	{
		int64_t *temp = (int64_t *)array.data();
		std::vector<int64_t> ret(temp, temp + array.size()/sizeof(int64_t));
		return ret;
	}
	
	std::vector<unsigned char> ArrayEncoder::Encode(std::vector<int64_t> array)
	{
		
 		// assert(array != nullptr);
		unsigned char *ptr = (unsigned char *)array.data();
		
		std::vector<unsigned char> newArray(ptr, ptr + array.size()*sizeof(int64_t));

		return newArray;
	}
	
	std::vector<unsigned char> ArrayEncoder::Encode(std::vector<uint64_t> array)
	{
		
 		// assert(array != nullptr);
		unsigned char *ptr = (unsigned char *)array.data();
		
		std::vector<unsigned char> newArray(ptr, ptr + array.size()*sizeof(uint64_t));

		return newArray;
	}

	std::vector<unsigned char> ArrayEncoder::Hash(std::vector<int64_t> array)
	{
 		// MyDebug::NonNull({array});
 		// assert(array->Length != 0);
		unsigned char *lengthInByte = int32ToByteArray(array.size());
		std::vector<unsigned char> lengthVec(lengthInByte, lengthInByte + sizeof(int32_t));
		
		std::vector<unsigned char> encoding = Encode(array);
		std::vector<unsigned char> hash1 = CryptoUtility::ComputeHash(encoding);		
		lengthVec.insert(lengthVec.end(), hash1.begin(), hash1.end());
		
		std::vector<unsigned char> hash2 = CryptoUtility::ComputeHash(lengthVec);
		
		return hash2;
	}
	
	std::vector<unsigned char> ArrayEncoder::Hash(std::vector<uint64_t> array)
	{
 		// MyDebug::NonNull({array});
 		// assert(array->Length != 0);
		unsigned char *lengthInByte = int32ToByteArray(array.size());
		std::vector<unsigned char> lengthVec(lengthInByte, lengthInByte + sizeof(int32_t));
		
		std::vector<unsigned char> encoding = Encode(array);
		std::vector<unsigned char> hash1 = CryptoUtility::ComputeHash(encoding);		
		lengthVec.insert(lengthVec.end(), hash1.begin(), hash1.end());
		
		std::vector<unsigned char> hash2 = CryptoUtility::ComputeHash(lengthVec);
		
		return hash2;
	}

	std::vector<int64_t> ArrayEncoder::Decodeint64_t(std::vector<unsigned char> &array)
	{
		int64_t *temp = (int64_t *)array.data();
		std::vector<int64_t> ret(temp, temp + array.size()/sizeof(int64_t));
		return ret;
	}

	std::vector<unsigned char> ArrayEncoder::EncodeInt64Array(std::vector<std::vector<int64_t> > array)
	{
		// array to vector of unsigned char
		// structure: ||# of layers||# of mul gates per layers 0||...||# of mul gates per layers t||data of layer 0||...||data of layer t||
		int arraySize = 0;
		for(int idx = 0; idx < array.size(); idx++)
		{
			arraySize += array[idx].size();
		}
		
		int64_t *data = new int64_t[1 + array.size() + arraySize];
		
		// ||# of layers||
		data[0] = array.size();
		
		int startPoint = 1 + array.size();
		
		int64_t *dataPtr = data + startPoint;
		
		for(int idx = 0; idx < array.size(); idx++)
		{
			// ||# of mul gates per layers idx||
			data[idx + 1] = array[idx].size();
			
			// ||data of layer idx||
			for(int kdx = 0; kdx < array[idx].size(); kdx++)
			{
				*dataPtr++ = array[idx][kdx];
			}
		}
		
		unsigned char *temp = (unsigned char *)data;
		
		std::vector<unsigned char> tempVec; 
		tempVec.insert(tempVec.end(), temp, temp + (1 + array.size() + arraySize)*sizeof(int64_t));
		return tempVec;
	}
	
	std::vector<unsigned char> ArrayEncoder::EncodeUInt64Array(std::vector<std::vector<uint64_t> > array)
	{
		// array to vector of unsigned char
		// structure: ||# of layers||# of mul gates per layers 0||...||# of mul gates per layers t||data of layer 0||...||data of layer t||
		int arraySize = 0;
		for(int idx = 0; idx < array.size(); idx++)
		{
			arraySize += array[idx].size();
		}
		
		uint64_t *data = new uint64_t[1 + array.size() + arraySize];
		
		// ||# of layers||
		data[0] = array.size();
		
		int startPoint = 1 + array.size();
		
		uint64_t *dataPtr = data + startPoint;
		
		for(int idx = 0; idx < array.size(); idx++)
		{
			// ||# of mul gates per layers idx||
			data[idx + 1] = array[idx].size();
			
			// ||data of layer idx||
			for(int kdx = 0; kdx < array[idx].size(); kdx++)
			{
				*dataPtr++ = array[idx][kdx];
			}
		}
		
		unsigned char *temp = (unsigned char *)data;
		
		std::vector<unsigned char> tempVec; 
		tempVec.insert(tempVec.end(), temp, temp + (1 + array.size() + arraySize)*sizeof(uint64_t));
		return tempVec;
	}

	std::vector<std::vector<int64_t> > ArrayEncoder::DecodeInt64Array(std::vector<unsigned char> array)
	{
 		// std::cout << "DecodeInt64Array: ";
		int64_t *data = (int64_t *)(array.data());
		int count = data[0];
 		// std::cout << data[0] << " ";
		std::vector<int> mulGatesPerLayer;
		
		for(int idx = 0; idx < count; idx++)
		{
			mulGatesPerLayer.push_back(data[idx + 1]);
 			// std::cout << data[idx + 1] << " ";
		}
		
		int startPoint = 1 + mulGatesPerLayer.size();;
		
		std::vector<std::vector<int64_t> > ret(0);
		
		for(int idx = 0; idx < count; idx++)
		{
			std::vector<int64_t> temp(mulGatesPerLayer[idx]);
			int64_t *dataPtr = data + startPoint;
			
			for(int kdx = 0; kdx < mulGatesPerLayer[idx]; kdx++)
			{
				temp[kdx] = *dataPtr++;
			}
			
			startPoint += mulGatesPerLayer[idx];
			ret.push_back(temp);
		}
		
		return ret;
	}
	
	std::vector<std::vector<uint64_t> > ArrayEncoder::DecodeUInt64Array(std::vector<unsigned char> array)
	{
 		// std::cout << "DecodeInt64Array: ";
		uint64_t *data = (uint64_t *)(array.data());
		int count = data[0];
 		// std::cout << data[0] << " ";
		std::vector<int> mulGatesPerLayer;
		
		for(int idx = 0; idx < count; idx++)
		{
			mulGatesPerLayer.push_back(data[idx + 1]);
 			// std::cout << data[idx + 1] << " ";
		}
		
		int startPoint = 1 + mulGatesPerLayer.size();;
		
		std::vector<std::vector<uint64_t> > ret(0);
		
		for(int idx = 0; idx < count; idx++)
		{
			std::vector<uint64_t> temp(mulGatesPerLayer[idx]);
			uint64_t *dataPtr = data + startPoint;
			
			for(int kdx = 0; kdx < mulGatesPerLayer[idx]; kdx++)
			{
				temp[kdx] = *dataPtr++;
			}
			
			startPoint += mulGatesPerLayer[idx];
			ret.push_back(temp);
		}
		
		return ret;
	}
	
	int ArrayEncoder::byteArrayToInt32(unsigned char *byteArray)
	{
		int val = 0;
		val |= (byteArray[0] << 24); 
		val |= (byteArray[1] << 16); 
		val |= (byteArray[2] << 8); 
		val |= (byteArray[3]); 
		return val;
	}
	
	unsigned char * ArrayEncoder::int32ToByteArray(unsigned int val)
	{
		unsigned char *byteArray = new unsigned char[4];

		byteArray[0] = val >> 24;
		byteArray[1] = val >> 16;
		byteArray[2] = val >> 8;
		byteArray[3] = val >> 0;
		
		return byteArray;
	}

	std::vector<unsigned char> ArrayEncoder::MacVectorToByteArray(std::vector<uint64_t> MACedEdges)
	{
		// 40 bit macs is equivalent to 5 bytes 
		std::vector<unsigned char> byteArray(5*MACedEdges.size());
		
		std::cout << "Converting Mac to byte array" << std::endl;
		
		int counter = 0;
		
		for(int idx = 0; idx < MACedEdges.size(); idx++)
		{
			unsigned char *data = (unsigned char *)(&MACedEdges[idx]);
			
			for(int kdx = 0; kdx < 5; kdx++)
			{
				byteArray[counter++] = data[kdx];
			}
		}
		
		return byteArray;
	}
	
	std::vector<unsigned char> ArrayEncoder::MacVectorToByteArray(uint64_t *MACedEdges, int size)
	{
		// 40 bit macs is equivalent to 5 bytes 
		std::vector<unsigned char> byteArray(5*size);
		
		std::cout << "Converting Mac to byte array" << std::endl;
		
		int counter = 0;
		
		for(int idx = 0; idx < size; idx++)
		{
			unsigned char *data = (unsigned char *)(&MACedEdges[idx]);
			
			for(int kdx = 0; kdx < 5; kdx++)
			{
				byteArray[counter] = data[kdx]; counter++;
			}
		}
		
		return byteArray;
	}
		
	std::vector<uint64_t> ArrayEncoder::ByteArrayToMacVector(std::vector<unsigned char> byteArray)
	{
		assert(byteArray.size() % 5 == 0);
		
		std::vector<uint64_t> MACedEdges(byteArray.size()/5);
		
		int counter = 0;
		for(int idx = 0; idx < MACedEdges.size(); idx++)
		{
			MACedEdges[idx] = 0;
			
			unsigned char *data = (unsigned char *)(&MACedEdges[idx]);
			for(int kdx = 0; kdx < 5; kdx++)
			{
				data[kdx] = byteArray[counter++];
			}
			data[5] = 0;
			data[6] = 0;
			data[7] = 0;
		}
		
		return MACedEdges;
	}
	
	void ArrayEncoder::ZZArray2ByteArray(unsigned char *dest, const ZZ *src, int nBytes, int size)
	{
		std::memset(dest, 0, nBytes*size);
		unsigned char *ptr = dest;
		for(int idx = 0; idx < size; idx++)
		{
			BytesFromZZ(ptr, src[idx], nBytes);
			ptr += nBytes;
		}
	}
		
	void ArrayEncoder::ZZArray2ByteArray(unsigned char *dest, const ZZ_p *src, int nBytes, int size)
	{
		std::memset(dest, 0, nBytes*size);
		unsigned char *ptr = dest;
		for(int idx = 0; idx < size; idx++)
		{
			BytesFromZZ(ptr, rep(src[idx]), nBytes);
			ptr += nBytes;
		}
	}
	
	void ArrayEncoder::ZZpX2ByteArray(unsigned char *dest, const ZZ_pX& src, int nBytes, int size)
	{
		std::memset(dest, 0, nBytes*size);
		unsigned char *ptr = dest;
		for(int idx = 0; idx < size; idx++)
		{
			BytesFromZZ(ptr, rep(src[idx]), nBytes);
			ptr += nBytes;
		}
	}
	
	void ArrayEncoder::ByteArray2ZZpX(ZZ_pX& dest, unsigned char *src, int nBytes, int size)
	{
		dest.SetLength(size);
		
		unsigned char *ptr = src;
		ZZ temp;
		
		for(int idx = 0; idx < size; idx++)
		{
			ZZFromBytes(temp, ptr, nBytes);
			ptr += nBytes;
			conv(dest[idx], temp);
		}
	}
	
	void ArrayEncoder::ByteArray2ZZArray(ZZ *dest, unsigned char *src, int nBytes, int size)
	{
		unsigned char *ptr = src;
		for(int idx = 0; idx < size; idx++)
		{
			ZZFromBytes(dest[idx], ptr, nBytes);
			ptr += nBytes;
		}
	}
		
	void ArrayEncoder::ByteArray2ZZArray(ZZ_p *dest, unsigned char *src, int nBytes, int size)
	{
		unsigned char *ptr = src;
		
		ZZ temp;
		
		for(int idx = 0; idx < size; idx++)
		{
			ZZFromBytes(temp, ptr, nBytes);
			ptr += nBytes;
			conv(dest[idx], temp);
		}
	}
	      
	DistributedMACGenerator::DistributedMACGenerator(uint64_t alpha, uint64_t beta) : alpha(alpha), beta(beta)
	{
		uint64_t val = 1;
		for(int idx = 0; idx < 50; idx++)
		{
			val = (val*alpha) % p;
			coeff.push_back(val);
		}
	}


	uint64_t DistributedMACGenerator::computeMAC(int party, std::vector<uint64_t> &data, uint64_t random)
	{
		uint64_t mac = 0;
		
		for(int idx = 0; idx < data.size(); idx++)
		{
			mac += data[idx]*coeff[idx];
		}
		
		if(party % 2)
		{
			mac += (beta + random);
		}
		else // Bob subtracts random from the mac
		{
			mac = mac - random;
		}
		
// 		return mac;
// 		return (uint64_t) (mac % Z64);
		return (mac % p);
	}
	
	SPDZ2kMAC::SPDZ2kMAC(uint64_t alpha, Communicator *com) : alpha(alpha), com(com){}
	
	uint64_t SPDZ2kMAC::computeMAC(int party, uint64_t data, uint64_t random)
	{
		if(party % 2 == 0)
		{
			return (alpha*data + random);
		}
		else
		{
			return (alpha*data - random);
		}
	}
	
	uint64_t SPDZ2kMAC::singleCheck(int party, uint64_t mac, uint64_t share)
	{
		uint64_t ourShare = mac - alpha*share;
		uint64_t theirShare;
		
		if(party % 2 == 0){
			com->SendEvaluationPartner((unsigned char *)(&ourShare), sizeof(uint64_t));
		}
		else
		{
			com->AwaitEvaluationPartner((unsigned char *)(&theirShare), sizeof(uint64_t));
		}
		
		// verify that the sum is 0
		return (ourShare + theirShare) == 0;
	}
	
// 	uint64_t SPDZ2kMAC::multipleCheck(int party, std::vector<uint64_t> mac, std::vector<uint64_t> share)
// 	{
// 		return 0;
// 	} 
	
	void TestUtility::PrintByteArray(const std::vector<unsigned char>& array, const std::string& str)
	{
		std::stringstream ss;
		ss << str << "\n";
		for(int idx = 0; idx < array.size(); idx++)
		{
			ss << (int)(array[idx]) << " ";
		}
		ss << "\n";
		
		std::cout << ss.str() << std::endl;
	}
	
	void TestUtility::PrintVector(const std::vector<int16_t>& input, const std::string& str){
		std::stringstream ss;
		ss << str << "\n";
		for(int idx = 0; idx < input.size(); idx++)
		{
			ss << (int)(input[idx]) << " ";
		}
		ss << "\n";
		
		std::cout << ss.str() << std::endl;
	}
	
	void TestUtility::PrintVector(const std::vector<uint16_t>& input, const std::string& str){
		std::stringstream ss;
		ss << str << "\n";
		for(int idx = 0; idx < input.size(); idx++)
		{
			ss << (int)(input[idx]) << " ";
		}
		ss << "\n";
		
		std::cout << ss.str() << std::endl;
	}
	
	void TestUtility::PrintVector(const std::vector<int64_t>& input, const std::string& str){
		std::stringstream ss;
		ss << str << "\n";
		for(int idx = 0; idx < input.size(); idx++)
		{
			ss << (int)(input[idx]) << " ";
		}
		ss << "\n";
		
		std::cout << ss.str() << std::endl;
	}
	
	void TestUtility::PrintVector(const std::vector<uint64_t>& input, const std::string& str){
		std::stringstream ss;
		ss << str << "\n";
		for(int idx = 0; idx < input.size(); idx++)
		{
			ss << (int)(input[idx]) << " ";
		}
		ss << "\n";
		
		std::cout << ss.str() << std::endl;
	}
	
	void TestUtility::PrintVector(const std::vector<std::vector<int64_t> >& input, const std::string& str, float scale){
		std::stringstream ss;
		ss << str << "\n";
		for(int idx = 0; idx < input.size(); idx++)
		{
			ss << idx << ": ";
			for(int kdx = 0; kdx < input[idx].size(); kdx++)
			{
				ss << (float)(input[idx][kdx]/scale) << " ";
			}
			ss << "\n";
		}
		
		std::cout << ss.str() << std::endl;
	}
	
	void TestUtility::PrintVector(const std::vector<std::vector<uint64_t> >& input, const std::string& str, float scale){
		std::stringstream ss;
		ss << str << "\n";
		for(int idx = 0; idx < input.size(); idx++)
		{
			ss << idx << ": ";
			for(int kdx = 0; kdx < input[idx].size(); kdx++)
			{
				ss << (float)(input[idx][kdx]/scale) << " ";
			}
			ss << "\n";
		}
		
		std::cout << ss.str() << std::endl;
	}
	
	void TestUtility::PrintMask(const std::vector<int64_t>& masks, const std::vector<int>& maskIndex, const std::vector<int>& itemsPerUser, const std::string& str){
		// Count number of required random values: some wires need 10 masks, some 1
		int numUsers = itemsPerUser.size();
		int numEdges = 0;
		for(int idx = 0; idx < itemsPerUser.size(); idx++)
		{
			numEdges += itemsPerUser[idx];
		}
		
		int numWires = 6*numUsers + 7*numEdges;		
		std::stringstream ss;
		ss << str << "\n";
		for(int idx = 0; idx < numWires-1; idx++)
		{
			ss << idx << ": ";
			for(int kdx = maskIndex[idx]; kdx < maskIndex[idx+1]; kdx++)
			{
				ss << (int64_t)(masks[kdx]) << " ";
			}
			ss << "\n";
		}
		
		for(int kdx = maskIndex[numWires-1]; kdx < masks.size(); kdx++)
		{
			ss << (int64_t)(masks[kdx]) << " ";
		}
		ss << "\n";
		
		std::cout << ss.str() << std::endl;
	}
	
	std::vector<std::vector<int64_t> > TestUtility::GenerateInput(const std::vector<int>& itemsPerUser, int inputLength)
	{
		unsigned seed =  std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
		std::uniform_int_distribution<int> distribution(-31, 31);
		std::uniform_int_distribution<int> distribution2(1, 5);
		
		std::vector<std::vector<int64_t> > input(inputLength);
		int count = 0;
		for(int udx = 0; udx < itemsPerUser.size(); udx++)
		{
			// Input for u
			for(int kdx = 0; kdx < 10; kdx++)
			{
			      int64_t rng = ((distribution(generator)) << 17);
			      
			      input[count].push_back(rng);
			}
			count++;
			// Input for items
			for(int idx = 0; idx < itemsPerUser[udx]; idx++)
			{
				// Input for v
				for(int kdx = 0; kdx < 10; kdx++)
				{
				      int64_t rng = ((distribution(generator)) << 17);
				      
				      input[count].push_back(rng);
				}
				count++;
				
				// Rating
				input[count].push_back(distribution2(generator)<<20);
				count++;
				
				// Is Real?
				input[count].push_back(1<<20);
				count++;
			}
		}
		
		return input;
	}
}
