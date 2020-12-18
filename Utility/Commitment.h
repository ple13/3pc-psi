#ifndef COMMITMENT_H__
#define COMMITMENT_H__

#pragma once

#include <string>
#include <vector>
#include <stdexcept>

#define SEED_LENGTH 64
#define COMMITMENT_LENGTH 32

namespace Utility
{
	/// <summary>
	/// Exception to throw when decommitment information does not match commitment.
	/// </summary>
	class InvalidDecommitmentException : public std::exception
	{
	public:
		InvalidDecommitmentException(const std::string& str);
	};

	/// <summary>
	/// Represents a commitment scheme.
	/// </summary>
	class ICommitment
	{
	public:
		virtual std::vector<unsigned char> Commit(std::vector<unsigned char> msg, std::vector<unsigned char>& seed) = 0;
		virtual bool Verification(std::vector<unsigned char> commitment, std::vector<unsigned char> msg, std::vector<unsigned char> seed) = 0;
	};

	class HashCommitment : public ICommitment
	{
	public:
		std::vector<unsigned char> CreateCommitment(std::vector<unsigned char> msg, std::vector<unsigned char> seed);

		std::vector<unsigned char> Commit(std::vector<unsigned char> msg, std::vector<unsigned char>& seed);

		bool Verification(std::vector<unsigned char> commitment, std::vector<unsigned char> msg, std::vector<unsigned char> seed);
	};
}

#endif
