
#ifndef ISECURE_RNG_H___
#define ISECURE_RNG_H___

#pragma once

#include <vector>
#include <stdexcept>
#include <random>
#include <ctime>
#include <chrono>

#include <openssl/aes.h>
#include <openssl/rand.h>
#include <openssl/conf.h>
#include <openssl/evp.h>
#include <openssl/err.h>

#include <cryptopp/rsa.h>
#include <cryptopp/aes.h>
#include "cryptopp/modes.h"
#include "cryptopp/filters.h"
#include <cryptopp/integer.h>
#include "cryptopp/hex.h"
#include <cryptopp/osrng.h>

#include <NTL/ZZ_p.h>

using namespace NTL;

namespace Utility
{
	class AESRNG// : public ISecureRNG
	{
	private:
		uint64_t counter;
		unsigned char * key;
		CryptoPP::ECB_Mode< CryptoPP::AES >::Encryption enc;
		// Return the cipher text
		unsigned char * Generate(unsigned char * plaintext, int plaintext_len, unsigned char * key);

	public:
		AESRNG();
		AESRNG(unsigned char *seed);
		
		~AESRNG()
		{
 			// delete [] key;
		}

		unsigned char * GetByteArray();
		unsigned __int128 GetByteArray(uint64_t counter);
// 		unsigned __int128 GetByteArray(unsigned char *data, int size);
		unsigned char * GetByteArray(unsigned char *data, int size);
		std::vector<uint16_t> GetUInt16Array(int length);
		std::vector<uint32_t> GetUInt32Array(int length);
		std::vector<uint64_t> GetUInt64Array(int length);
		std::vector<uint64_t> GetUInt64Array(__uint128_t key, int length);
		std::vector<unsigned __int128> GetUInt128Array(int length);
		void ComputeFx(std::vector<ZZ_p>& dataFx, const std::vector<uint64_t>& data);
		void ComputeRingFx(std::vector<__uint128_t>& res, const std::vector<uint64_t>& data);
		void ComputeRingFx(std::vector<__uint128_t>& res, const std::vector<__uint128_t>& data);
		void ComputeFx(std::vector<ZZ_p>& dataFx, const std::vector<ZZ_p>& data, int nBytes);
		void GetRandomFx(std::vector<ZZ_p>& res, int nBytes, int length);
		void GetRandomRingFx(std::vector<__uint128_t>& res, int nBytes, int length);
		std::vector<int64_t>  GetMaskArray(int length);
		
		unsigned char * AES(std::string plain);
	};
}

#endif
