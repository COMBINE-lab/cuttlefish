#pragma once
#include "parallel-queues-common.h"
#include <string>
#include <sstream>
#include <set>
#include <mutex>
#include <atomic>
#include <vector>
#include <iostream>

namespace refresh {

	class queue_monitor
	{
		std::ostream& ostr;
		bool single_line;
		bool reporting;
		std::mutex mtx;

		std::vector<IQueuePrinter*> queue_printers;
	public:
		queue_monitor(
			std::ostream& ostr,
			bool single_line = false,
			bool _reporting = true);
		void enable_reporting() { reporting = true; }
		void disable_reporting() { reporting = false; }
		void report();
		void force_report();
		void print_summary();
		void register_queue(IQueuePrinter* queue_printer);
	};

	class basic_queue_observer : public IQueueObserver, public IQueuePrinter {
		queue_monitor& monitor;
		std::string name;
		size_t max_size = std::numeric_limits<size_t>::max();
		size_t cur_size{};
		bool completed = false;

		std::chrono::nanoseconds tot_pop_wait_time{};
		std::chrono::nanoseconds tot_push_wait_time{};

	public:
		basic_queue_observer(
			queue_monitor& monitor) :
			monitor(monitor)
		{
			monitor.register_queue(this);
		}

		void set_queue_params(const std::string& name, size_t max_size) override
		{
			this->name = name;
			this->max_size = max_size;
		}
		void notify_pushed() override
		{
			++cur_size;
			monitor.report();
		}
		void notify_popped() override
		{
			--cur_size;
			monitor.report();
		}

		void notify_completed() override
		{
			completed = true;
		}

		void print(std::ostream& oss, const std::string& sep) const override
		{
			//dont print empty completed queues
			if (completed && cur_size == 0)
				return;

			oss << name;
			if (completed)
				oss << " (c)";
			oss << ": ";
			if (max_size == std::numeric_limits<size_t>::max()) //skip max_size for unlimited queues
				oss << cur_size;
			else
				oss << cur_size << " " << max_size;

			oss << sep;
		}
		void print_summary(std::ostream& oss) const override
		{
			oss << name << ":\n";
			oss << "\twait on pop tot time: " << std::chrono::duration<double>(tot_pop_wait_time).count() << "s\n";
			oss << "\twait on push tot time: " << std::chrono::duration<double>(tot_push_wait_time).count() << "s\n";
		}

		void notify_wait_on_push_time(std::chrono::nanoseconds time) override
		{
			tot_push_wait_time += time;
		}
		void notify_wait_on_pop_time(std::chrono::nanoseconds time) override
		{
			tot_pop_wait_time += time;
		}
	};

} //namespace refresh