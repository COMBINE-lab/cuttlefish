#ifndef _KEYBOARD_LISTENER_H
#define _KEYBOARD_LISTENER_H

#include <iostream>
#include <thread>
#include <mutex>
#include <functional>
#include <map>
#include <atomic>
#ifdef _WIN32
#include <windows.h>
#endif

class key_listener {
    std::thread th;
    std::atomic_bool stop = false;
    using callback_t = std::function<void(const std::string& /*command*/)>;
    std::map<std::string, callback_t> observers;
    std::mutex mtx;

    key_listener() :
    th([this]{
        std::string name, command;
        while(!stop) {
            std::cin >> name;
            auto it = observers.find(name);
            if(it == observers.end()) {
                std::cerr << "There is no '" << name << "' observer\n";
            } else {
                std::cin >> command;
                it->second(command);
            }
        }
    })
    {
    }

    public:
    static key_listener& Inst() {
        static key_listener inst;
        return inst;
    }
    void register_observer(const std::string& name, callback_t ob) {
        std::lock_guard<std::mutex> lck(mtx);
        if (observers.find(name) != observers.end()) {
            std::cerr << "Observer of name " << name << " already exists, skipping\n";
            return;
        }
        std::cerr << "Registering observer " << name << "\n";
        observers.emplace(name, ob);
    }
    void unregister_observer(const std::string& name) {
        std::lock_guard<std::mutex> lck(mtx);
        auto it = observers.find(name);
        if (it == observers.end()) {
            std::cerr << "Observer " << name << " was not registered, cannot be unregitered\n";
            return;
        }
        observers.erase(it);
    }

    ~key_listener() {
        stop = true;
        auto hande = th.native_handle();
#ifdef _WIN32
        TerminateThread(hande, 0);
#else
        pthread_cancel(hande);
#endif
        //std::cerr << "To terminate you must enter anything...\n";
        th.join();
    }
};

#endif // _KEYBOARD_LISTENER_H
