#include <cassert>
#include <stdio.h>
#include <string>
#include <sstream>
#include <iostream>
#include "Channels.h"


namespace Utility
{
	namespace Properties
	{

		void OneWayCommunicationChannel::SendMessage(const std::vector<unsigned char>& msg)
		{
			assert(msg.size() > 0);
			while (message.size() > 0)
			{
			}
			
			this->message = msg;
		}

		void OneWayCommunicationChannel::SendMessage(std::vector<unsigned char>&& msg)
		{
			assert(msg.size() > 0);
			while (message.size() > 0)
			{
			}
			
			this->message = std::move(msg);
		}
		
		std::vector<unsigned char> OneWayCommunicationChannel::AwaitMessage()
		{
			while (message.empty())
			{
			}
			
			auto temp = message;
			message.resize(0);
			return temp;
		}

		Sender::Sender(OneWayCommunicationChannel *channel) : channel(channel)
		{
		}

		void Sender::SendMessage(const std::vector<unsigned char>& msg)
		{
			channel->SendMessage(msg);
		}
		
		void Sender::SendMessage(std::vector<unsigned char>&& msg)
		{
			channel->SendMessage(std::move(msg));
		}

		Receiver::Receiver(OneWayCommunicationChannel *channel) : channel(channel)
		{
		}

		std::vector<unsigned char> Receiver::AwaitMessage()
		{
			return channel->AwaitMessage();
		}
	}
}
