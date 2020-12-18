#pragma once

#include <vector>

namespace Utility
{
	class ICommunicator
	{
	public:
		virtual std::vector<unsigned char>  AwaitAlice() = 0;
		virtual std::vector<unsigned char> AwaitBob() = 0;
		virtual std::vector<unsigned char>  AwaitEvaluationPartner() = 0;
		void AwaitFromPlayers(std::vector<unsigned char>  &m1, std::vector<unsigned char> &m2, std::vector<unsigned char> &m3);
		virtual std::vector<unsigned char> AwaitMessage(int player) = 0;
		virtual std::vector<unsigned char> AwaitNonPartner() = 0;
		virtual std::vector<unsigned char> AwaitVerificationPartner() = 0;
		
		virtual void SendAlice(const std::vector<unsigned char>& msg) = 0;
		virtual void SendAlice(std::vector<unsigned char>&& msg) = 0;
		
		virtual void SendBob(const std::vector<unsigned char>& msg) = 0;
		virtual void SendBob(std::vector<unsigned char>&& msg) = 0;
		
		virtual void SendEvaluationPartner(const std::vector<unsigned char>& msg) = 0;
		virtual void SendEvaluationPartner(std::vector<unsigned char>&& msg) = 0;
		
		virtual void SendMessage(int player, const std::vector<unsigned char>& msg) = 0;
		virtual void SendMessage(int player, std::vector<unsigned char>&& msg) = 0;
		
		virtual void SendNonPartner(const std::vector<unsigned char>& msg) = 0;
		virtual void SendNonPartner(std::vector<unsigned char>&& msg) = 0;
		
		virtual void SendToAll(const std::vector<unsigned char>& m) = 0;
		virtual void SendToPlayers(const std::vector<unsigned char>& m1, const std::vector<unsigned char>& m2, const std::vector<unsigned char>& m3) = 0;
		
		virtual void SendVerificationPartner(const std::vector<unsigned char>& msg) = 0;
		virtual void SendVerificationPartner(std::vector<unsigned char>&& msg) = 0;
		
		virtual void Transfer(const std::vector<unsigned char>& msg) = 0;
	};
}
