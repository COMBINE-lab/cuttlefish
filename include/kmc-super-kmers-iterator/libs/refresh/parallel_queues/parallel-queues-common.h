#ifndef _REFRESH_PARALLEL_QUEUES_COMMON
#define _REFRESH_PARALLEL_QUEUES_COMMON

#include <ostream>
#include <chrono>
#include <string>

//uncomment for queue monitoring
//#define REFRESH_PROFILE_QUEUES_DETAILED
//#define REFRESH_PROFILE_QUEUES
//#define REFRESH_USE_KEYBOARD_LISTENER_FOR_PARALLEL_QUEUES

namespace refresh {
	class IQueueObserver {
	public:
		virtual void set_queue_params(const std::string& name, size_t max_size) = 0;
		virtual void notify_pushed() = 0;
		virtual void notify_popped() = 0;
		virtual void notify_completed() = 0;
		virtual void notify_wait_on_push_time(std::chrono::nanoseconds time) = 0;
		virtual void notify_wait_on_pop_time(std::chrono::nanoseconds time) = 0;
		virtual ~IQueueObserver() = default;
	};

	class IQueuePrinter {
	public:
		virtual void print(std::ostream& oss, const std::string& sep) const = 0;
		virtual void print_summary(std::ostream& oss) const = 0;
		virtual ~IQueuePrinter() = default;
	};

#ifdef REFRESH_PROFILE_QUEUES
	refresh::IQueueObserver* refresh_profile_queues_create_observer();
	std::string refresh_profile_queues_create_queue_name();
#endif 

} //namespace refresh

#endif // _REFRESH_PARALLEL_QUEUES_COMMON
