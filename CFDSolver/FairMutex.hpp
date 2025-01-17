#pragma once

#include <thread>
#include <mutex>
#include <optional>
#include <queue>
#include <condition_variable>

class FairMutex {
public:
	void lock() {
		std::unique_lock<std::mutex> lock(mutex);
		std::condition_variable cv;
		wait_queue.push(&cv);
		cv.wait(lock, [&] { return wait_queue.front() == &cv && !locked; });
		locked = true;
		wait_queue.pop();
	}

	void unlock() {
		std::lock_guard<std::mutex> lock(mutex);
		locked = false;
		if (!wait_queue.empty())
			wait_queue.front()->notify_one();
	}

private:
	std::mutex mutex;
	std::queue<std::condition_variable*> wait_queue;
	bool locked = false;
};