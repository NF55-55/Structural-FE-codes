#ifndef TIMER_HIGH_RES_H
#define TIMER_HIGH_RES_H

#include <iostream>
#include <chrono>
#include <string_view>

class Timer
{
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start {};
    std::chrono::time_point<std::chrono::high_resolution_clock> m_stop {};
    bool m_active {true};

public:
    Timer() = default;

    void start(){ m_start = std::chrono::high_resolution_clock::now(); }
    void stop(){ m_stop = std::chrono::high_resolution_clock::now(); }

    void setOn(){ m_active = true; }
    void setOff(){ m_active = false; }

    auto getDurationMicro() const{
        std::chrono::duration<double,std::milli> d {(m_stop-m_start)*1000};
        return d.count();
    }

    void printDurationMicro(const int refNum
        , const std::string_view units="micro s", const int divisor=1) const
    {
        if (m_active){
            auto durationMicro {getDurationMicro()};
            std::cout << "Time elapsed(" << refNum << "): " 
                << durationMicro/divisor << ' ' << units << '\n';
        }
    }
};

#endif