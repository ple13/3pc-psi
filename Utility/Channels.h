#ifndef CHANNEL_H__
#define CHANNEL_H__

#pragma once

#include <vector>
#include <thread>
#include <mutex>
 

namespace Utility
{
	namespace Properties
	{
		class IInputChannel
		{
		public:
			virtual void SendMessage(const std::vector<unsigned char>& msg) = 0;
			virtual void SendMessage(std::vector<unsigned char>&& msg) = 0;
		};

		class IOutputChannel
		{
		public:
			virtual std::vector<unsigned char> AwaitMessage() = 0;
		};

		class OneWayCommunicationChannel
		{
		private:
			std::vector<unsigned char> message;
			std::mutex mutex;
			//public volatile object readLock = new object();
			//public volatile object writeLock = new object();

		public:
			void SendMessage(const std::vector<unsigned char>& msg);
			void SendMessage(std::vector<unsigned char>&& msg);

			std::vector<unsigned char> AwaitMessage();
		};

		class Sender : public IInputChannel
		{
		private:
			OneWayCommunicationChannel *const channel;

		public:
			virtual ~Sender()
			{
				delete channel;
			}

			Sender(OneWayCommunicationChannel *channel);

			void SendMessage(const std::vector<unsigned char>& msg);
			void SendMessage(std::vector<unsigned char>&& msg);
		};

		class Receiver : public IOutputChannel
		{
		private:
			OneWayCommunicationChannel *const channel;

		public:
			virtual ~Receiver()
			{
				delete channel;
			}

			Receiver(OneWayCommunicationChannel *channel);

			std::vector<unsigned char> AwaitMessage();
		};

		class ChannelBuilder
		{
		public:
			static void CreateChannel(IInputChannel * &sender, IOutputChannel * &receiver)
			{
				auto channel = new OneWayCommunicationChannel();
				sender = new Sender(channel);
				receiver = new Receiver(channel);
			}
		};
	}
}

#endif
