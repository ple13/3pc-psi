#include <cassert>
#include <iostream>
#include <openssl/rand.h>
#include <openssl/aes.h>
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
#include "CryptoUtility.h"
#include "ISecureRNG.h"
#include "Timer.h"

#include <emp-tool/utils/block.h>
#include <emp-tool/utils/prp.h>
#include <emp-tool/utils/prg.h>

#if !defined (__AES__)
    #error "AES-NI instructions not enabled"
#endif

namespace Utility
{
	AESRNG::AESRNG(){}
	
	/// Initialize an AES encryptor/decryptor
	AESRNG::AESRNG(unsigned char * seed)
	{
		key = new unsigned char[32];
		for(int idx = 0; idx < 32; idx++)
		{
			key[idx] = seed[idx];
		}
		counter = 0;
		
		CryptoPP::SecByteBlock key(this->key, 16);
		
		enc.SetKey( key, key.size() );
	}
	
	unsigned char * AESRNG::AES(std::string plain)
	{

		std::string cipher, encoded;

		try
		{
		    // The StreamTransformationFilter adds padding
		    //  as required. ECB and CBC Mode must be padded
		    //  to the block size of the cipher.
		    CryptoPP::StringSource ss1( plain, true, 
			new CryptoPP::StreamTransformationFilter( enc,
			    new CryptoPP::StringSink( cipher )
			) // StreamTransformationFilter      
		    ); // StringSource
		}
		catch( CryptoPP::Exception& e )
		{
		    cerr << e.what() << endl;
		    exit(1);
		}
		
// 		// Pretty print cipher text
		CryptoPP::StringSource ss2( cipher, true,
		    new CryptoPP::HexEncoder(
			new CryptoPP::StringSink( encoded )
		    ) // HexEncoder
		); // StringSource
		
		
		return (unsigned char *)(encoded.c_str());
// 		cout << "cipher text: " << encoded << endl;
	}

	unsigned char * AESRNG::GetByteArray()
	{
		counter++;
		unsigned char *temp = (unsigned char *)(&counter);
		
		return Generate(temp, sizeof(counter), key);
	}
	
	unsigned __int128 AESRNG::GetByteArray(uint64_t counter)
	{
		unsigned char *temp = (unsigned char *)(&counter);
		
		unsigned __int128 *val = (unsigned __int128 *)Generate(temp, sizeof(counter), key);
		
		return *val;
	}
	
	unsigned char * AESRNG::GetByteArray(unsigned char *data, int size)
	{		
		return Generate(data, size, key);
	}
	
	std::vector<uint16_t> AESRNG::GetUInt16Array(int length)
	{
		std::vector<uint16_t> ret(length);
		
		for(int idx = 0; idx < length/8; idx++)
		{
			uint16_t *prng = (uint16_t *)GetByteArray();
			
			ret[8*idx] = prng[0];
			ret[8*idx + 1] = prng[1];
			ret[8*idx + 2] = prng[2];
			ret[8*idx + 3] = prng[3];
			ret[8*idx + 4] = prng[4];
			ret[8*idx + 5] = prng[5];
			ret[8*idx + 6] = prng[6];
			ret[8*idx + 7] = prng[7];
			
			delete [] prng;
		}
		
		if(length - 8*(length/8) > 0)
		{
		    uint16_t *prng = (uint16_t *)GetByteArray();

		    for(int idx = 0; idx < length - 8*(length/8); idx++)
		    {
			    ret[8*(length/8) + idx] = prng[idx];
		    }
		    delete [] prng;
		}
		
		return ret;
	}
	std::vector<uint32_t> AESRNG::GetUInt32Array(int length)
	{
		std::vector<uint32_t> ret(length);
		
		for(int idx = 0; idx < length/4; idx++)
		{
			uint32_t *prng = (uint32_t *)GetByteArray();
			
			ret[4*idx] = prng[0];
			ret[4*idx + 1] = prng[1];
			ret[4*idx + 2] = prng[2];
			ret[4*idx + 3] = prng[3];
			
			delete [] prng;
		}
		
		if(length - 4*(length/4) > 0)
		{
		    uint32_t *prng = (uint32_t *)GetByteArray();

		    for(int idx = 0; idx < length - 4*(length/4); idx++)
		    {
			    ret[4*(length/4) + idx] = prng[idx];
		    }
		    delete [] prng;
		}
		
		return ret;
	}
	
	std::vector<uint64_t> AESRNG::GetUInt64Array(int length)
	{
		std::vector<uint64_t> ret(length);
		
		for(int idx = 0; idx < length/2; idx++)
		{
			uint64_t *prng = (uint64_t *)GetByteArray();
			ret[2*idx] = prng[0];
			ret[2*idx + 1] = prng[1];
			
			delete [] prng;
		}
		
		if(length - 2*(length/2) > 0)
		{
		      uint64_t *prng = (uint64_t *)GetByteArray();
		      
		      for(int idx = 0; idx < length - 2*(length/2); idx++)
		      {
			      ret[2*(length/2) + idx] = prng[idx];
		      }
		      delete [] prng;
		}
		
		return ret;
	}
	
	std::vector<unsigned __int128> AESRNG::GetUInt128Array(int length)
	{
		std::vector<unsigned __int128> ret(length);
#if !defined (__AES__)
		for(int idx = 0; idx < length; idx++)
		{
			unsigned char *prng = GetByteArray();
			ret[idx] = *((unsigned int *)prng);
			
			delete [] prng;
		}
#else
		emp::block *data = new emp::block[length];
		emp::PRG prg(emp::fix_key);
		prg.random_block(data, length);
		
		uint64_t *temp;
		for(int idx = 0; idx < length; idx++)
		{
			temp = (uint64_t *)(data + idx);
			ret[idx] = temp[1];
			ret[idx] = (ret[idx] >> 64) + temp[0];
		}
		delete [] data;
#endif
		return ret;
	}
	
	std::vector<uint64_t> AESRNG::GetUInt64Array(__uint128_t _key, int length)
	{
		std::vector<unsigned uint64_t> ret(length);
#if !defined (__AES__)
		for(int idx = 0; idx < length; idx++)
		{
			unsigned char *prng = GetByteArray();
			ret[idx] = *((unsigned int *)prng);
			
			delete [] prng;
		}
#else
		emp::block *data = new emp::block[length/2];
		emp::PRG prg((const char *)(&_key));
		prg.random_block(data, length/2);
		
		uint64_t *temp = (uint64_t *)data;
		
		for(int idx = 0; idx < length; idx++)
		{
			ret[idx] = temp[idx];
		}
		delete [] data;
#endif
		return ret;
	}
	
	void AESRNG::ComputeFx(std::vector<ZZ_p>& res, const std::vector<uint64_t>& data)
	{
		res.resize(data.size());
#if !defined (__AES__)
		ZZ temp;
		for(int idx = 0; idx < data.size(); idx++)
		{
			unsigned char *val = GetByteArray((unsigned char *)&data[idx], sizeof(uint64_t)); 
			ZZFromBytes(temp, val, 16);
			conv(res[idx], temp);
			delete val;
		}
#else
		emp::PRP prp(emp::fix_key);
		emp::block *temp = new emp::block[data.size()];
		
		for(int idx = 0; idx < data.size(); idx++)
		{
			temp[idx] = emp::makeBlock(0, data[idx]);
		}
		
		prp.permute_block(temp, data.size());
		
		for(int idx = 0; idx < data.size(); idx++)
		{
			conv(res[idx], ZZFromBytes((unsigned char *)(temp + idx), 16));
		}
		
		delete [] temp;
#endif
	}
	
	void AESRNG::ComputeRingFx(std::vector<__uint128_t>& res, const std::vector<uint64_t>& data)
	{
		__uint128_t modulo = 1; modulo = (modulo << 80) - 65;
		__uint128_t mm = 1; mm = (mm << 80) - 1;
		
		res.resize(data.size());
#if !defined (__AES__)
// 		ZZ temp;
// 		for(int idx = 0; idx < data.size(); idx++)
// 		{
// 			unsigned char *val = GetByteArray((unsigned char *)&data[idx], sizeof(uint64_t)); 
// 			ZZFromBytes(temp, val, 16);
// 			conv(res[idx], temp);
// 			delete val;
// 		}
#else
		emp::PRP prp(emp::fix_key);
		emp::block *temp = new emp::block[data.size()];
		
		for(int idx = 0; idx < data.size(); idx++)
		{
			temp[idx] = emp::makeBlock(0, data[idx]);
		}
		
		prp.permute_block(temp, data.size());
// 		t.Tick("permute_block");
		for(int idx = 0; idx < data.size(); idx++)
		{
			res[idx] = *((__uint128_t *)(temp + idx));
			res[idx] = (res[idx] >> 80) + ((res[idx] >> 80) << 6) + (res[idx] & mm);
			while(res[idx] >= modulo) res[idx] -= modulo;
		}
		
		delete [] temp;
#endif
	}
	
	void AESRNG::ComputeRingFx(std::vector<__uint128_t>& res, const std::vector<__uint128_t>& data)
	{
		__uint128_t modulo = 1; modulo = (modulo << 80) - 65;
		__uint128_t mm = 1; mm = (mm << 80) - 1;
		
		res.resize(data.size());
#if !defined (__AES__)
// 		ZZ temp;
// 		for(int idx = 0; idx < data.size(); idx++)
// 		{
// 			unsigned char *val = GetByteArray((unsigned char *)&data[idx], sizeof(uint64_t)); 
// 			ZZFromBytes(temp, val, 16);
// 			conv(res[idx], temp);
// 			delete val;
// 		}
#else
		emp::PRP prp(emp::fix_key);
		emp::block *temp = new emp::block[data.size()];
		
		for(int idx = 0; idx < data.size(); idx++)
		{
			temp[idx] = emp::makeBlock(0, data[idx]);
		}
		
		prp.permute_block(temp, data.size());
// 		t.Tick("permute_block");
		for(int idx = 0; idx < data.size(); idx++)
		{
			res[idx] = *((__uint128_t *)(temp + idx));
			res[idx] = (res[idx] >> 80) + ((res[idx] >> 80) << 6) + (res[idx] & mm);
			while(res[idx] >= modulo) res[idx] -= modulo;
		}
		
		delete [] temp;
#endif
	}
	
	void AESRNG::ComputeFx(std::vector<ZZ_p>& res, const std::vector<ZZ_p>& data, int nBytes)
	{
		res.resize(data.size());
		
#if !defined (__AES__)
		unsigned char *bytes = new unsigned char[nBytes*data.size()];
		unsigned char *ptr = bytes;
		
		ArrayEncoder::ZZArray2ByteArray(bytes, data.data(), nBytes, data.size());
		ZZ temp;
		for(int idx = 0; idx < data.size(); idx++)
		{
			unsigned char *val = GetByteArray(ptr, nBytes);
			conv(res[idx], ZZFromBytes(val, nBytes));
			ptr += nBytes;
		}
		
		ptr = 0;
		delete [] bytes;
#else
		unsigned char *bytes = new unsigned char[16*data.size()];
		
		ArrayEncoder::ZZArray2ByteArray(bytes, data.data(), 16, data.size());
		
		emp::block *temp = (emp::block *)bytes;
		
		emp::PRP prp(emp::fix_key);
		prp.permute_block(temp, data.size());
		
		for(int idx = 0; idx < data.size(); idx++)
		{
			conv(res[idx], ZZFromBytes((unsigned char *)(temp + idx), nBytes));
		}
		
		temp = 0;
		
		delete [] bytes;
#endif
	}
	
	void AESRNG::GetRandomFx(std::vector<ZZ_p>& res, int nBytes, int length)
	{
		res.resize(length);
#if !defined (__AES__)
		ZZ temp;
                for(int idx = 0; idx < length; idx++)
                {
                        unsigned char *val = GetByteArray((unsigned char *)&idx, sizeof(int));
                        conv(res[idx], ZZFromBytes(val, nBytes));
                }
#else
		Timer t;
		emp::block *data = new emp::block[length];
		emp::PRG prg(emp::fix_key);
		prg.random_block(data, length);
// 		t.Tick("Generating random numbers");
		for(int idx = 0; idx < length; idx++)
		{
			conv(res[idx], ZZFromBytes((unsigned char *)(data + idx), nBytes));
		}
// 		t.Tick("NTL converstion");
		delete [] data;
#endif
	}
	
	void AESRNG::GetRandomRingFx(std::vector<__uint128_t>& res, int nBytes, int length)
	{
		__uint128_t modulo = 1; modulo = (modulo << 80) - 65;
		res.resize(length);
#if !defined (__AES__)
// 		ZZ temp;
//                 for(int idx = 0; idx < length; idx++)
//                 {
//                         unsigned char *val = GetByteArray((unsigned char *)&idx, sizeof(int));
//                         conv(res[idx], ZZFromBytes(val, nBytes));
//                 }
#else
		Timer t;
		emp::block *data = new emp::block[length];
		emp::PRG prg(emp::fix_key);
		prg.random_block(data, length);
// 		t.Tick("Generating random numbers");
		for(int idx = 0; idx < length; idx++)
		{
			res[idx] = *((__uint128_t *)(data + idx)) >> 48;
		}
// 		t.Tick("Ring converstion");
		delete [] data;
#endif
	}
	
	std::vector<int64_t> AESRNG::GetMaskArray(int length)
	{
		std::vector<int64_t> ret(length);
		
		for(int idx = 0; idx < length/2; idx++)
		{
			int64_t *prng = (int64_t *)GetByteArray();
			
			ret[2*idx] = prng[0];
			ret[2*idx + 1] = prng[1];
			
			delete [] prng;
		}
		
		if(length - 2*(length/2) > 0)
		{
		      int64_t *prng = (int64_t *)GetByteArray();
		      
		      for(int idx = 0; idx < length - 2*(length/2); idx++)
		      {
			      ret[2*(length/2) + idx] = prng[idx];
		      }
		      delete [] prng;
		}
		
		return ret;
	}

	unsigned char * AESRNG::Generate(unsigned char * plaintext, int plaintext_len, unsigned char * key)
	{
		unsigned char *ciphertext = new unsigned char[16];
		int ciphertext_len;
		
		EVP_CIPHER_CTX *ctx;

		int len;

		/* Create and initialise the context */
		if(!(ctx = EVP_CIPHER_CTX_new())){
// 			handleErrors();
		}

		/* Initialise the encryption operation. IMPORTANT - ensure you use a key
		* In this example we are using 256 bit AES (i.e. a 256 bit key). 
		*/
		if(1 != EVP_EncryptInit_ex(ctx, EVP_aes_256_ecb(), NULL, key, NULL))
		{
// 			handleErrors();
		}

		/* Provide the message to be encrypted, and obtain the encrypted output.
		* EVP_EncryptUpdate can be called multiple times if necessary
		*/
		if(1 != EVP_EncryptUpdate(ctx, ciphertext, &len, plaintext, plaintext_len)) 
		{
// 			handleErrors();
		}
		
		ciphertext_len = len;

		/* Finalise the encryption. Further ciphertext bytes may be written at
		* this stage.
		*/
		if(1 != EVP_EncryptFinal_ex(ctx, ciphertext + len, &len))
		{
// 			handleErrors();
		}
		
		ciphertext_len += len;

		/* Clean up */
		EVP_CIPHER_CTX_free(ctx);
		
		return ciphertext;
	}
}
