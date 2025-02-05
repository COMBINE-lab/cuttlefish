#include "parallel-queues-monitor.h"
#ifdef REFRESH_USE_KEYBOARD_LISTENER_FOR_PARALLEL_QUEUES
#include "keyboard_listener.h"
#endif
#include <memory>

namespace refresh {
	queue_monitor::queue_monitor(
		std::ostream& ostr,
		bool single_line,
		bool _reporting)
		:
		ostr(ostr),
		single_line(single_line),
		reporting(_reporting)
	{
	}

	void queue_monitor::force_report() {
		std::lock_guard <std::mutex> lck(mtx);

		std::string sep = single_line ? "\t" : "\n";
		std::string rep;

		if (!single_line)
			rep = "**********\n";

		rep += std::to_string(std::time(nullptr)) + sep;

		std::ostringstream oss;
		for (auto x : queue_printers)
		{
			x->print(oss, sep);
		}

		ostr << oss.str();
		if (!rep.empty())
			rep.back() = '\n';

		if (!single_line)
			rep += "**********\n";

		ostr << rep;
	}
	void queue_monitor::report()
	{
		if (!reporting)
			return;

		force_report();

	}

	void queue_monitor::print_summary()
	{
		for (auto x : queue_printers)
			x->print_summary(ostr);
	}

	void queue_monitor::register_queue(IQueuePrinter* queue_printer)
	{
		std::lock_guard<std::mutex> lck(mtx);
		queue_printers.push_back(queue_printer);
	}

#ifdef REFRESH_PROFILE_QUEUES
	class refresh_profile_queues_global {
		queue_monitor monitor;
		std::vector<std::unique_ptr<basic_queue_observer>> observers;
		std::mutex mtx;
		refresh_profile_queues_global() :
#ifdef REFRESH_PROFILE_QUEUES_DETAILED
			monitor(std::cerr, true, true)
#else
			monitor(std::cerr, true, false)
#endif // REFRESH_PROFILE_QUEUES_DETAILED
		{
#ifdef REFRESH_USE_KEYBOARD_LISTENER_FOR_PARALLEL_QUEUES
			key_listener::Inst().register_observer("queues", [this](const std::string& command) {
				if (command == "report") {
					monitor.force_report();
				} else if (command == "enable") {
					monitor.enable_reporting();
				} else if (command == "disable") {
					monitor.disable_reporting();
				}
			});
#endif
		}


	public:
		static refresh_profile_queues_global& Inst() {
			static refresh_profile_queues_global inst;
			return inst;
		}

		IQueueObserver* create_observer()
		{
			std::lock_guard<std::mutex> lck(mtx);
			observers.push_back(std::make_unique<basic_queue_observer>(monitor));
			return observers.back().get();
		}
		~refresh_profile_queues_global() {
			monitor.print_summary();
#ifdef REFRESH_USE_KEYBOARD_LISTENER_FOR_PARALLEL_QUEUES
			key_listener::Inst().unregister_observer("queues");
#endif
		}
	};

	refresh::IQueueObserver* refresh_profile_queues_create_observer() {
		return refresh_profile_queues_global::Inst().create_observer();
	}
	std::string refresh_profile_queues_create_queue_name() {
		static std::atomic<size_t> n_calls{};
		return std::string("unnamed queue ") + std::to_string(n_calls++);
	}
#endif // REFRESH_PROFILE_QUEUES
} //namespace refresh


