#include "Commitment.h"
#include "CryptoUtility.h"
#include <openssl/sha.h>
#include <cstring>

namespace Utility
{

	InvalidDecommitmentException::InvalidDecommitmentException(const std::string& str) : std::exception(/*str*/)
	{
	}

	std::vector<unsigned char> HashCommitment::CreateCommitment(std::vector<unsigned char> msg, std::vector<unsigned char> seed)
	{
		msg.insert(msg.end(), seed.begin(), seed.end());
		
		return CryptoUtility::ComputeHash(msg);
	}

	std::vector<unsigned char> HashCommitment::Commit(std::vector<unsigned char> msg, std::vector<unsigned char>& seed)
	{
		seed = CryptoUtility::SampleByteArray(SEED_LENGTH);
		return CreateCommitment(msg, seed);
	}

	bool HashCommitment::Verification(std::vector<unsigned char> commitment, std::vector<unsigned char> msg, std::vector<unsigned char> seed)
	{
		return (commitment == CreateCommitment(msg, seed));
	}
}
