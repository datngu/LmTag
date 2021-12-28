#ifndef TIME_LOG
#define TIME_LOG
#include <ctime>
#include <chrono>

using namespace std::chrono;
class TimeLog
{
private:
    high_resolution_clock::time_point begin;
    system_clock::time_point sys_begin;
    high_resolution_clock::time_point end;
    system_clock::time_point sys_end;
    time_t tt;
public:

    TimeLog(){
        begin = high_resolution_clock::now();
        sys_begin = system_clock::now();
    };
    /*
    void begin_clock(){
        begin = high_resolution_clock::now();
        sys_begin = system_clock::now();
    };
    void end_clock(){
        end = high_resolution_clock::now();
        sys_end = system_clock::now();
    };
    */
    void print_clock(){
        end = high_resolution_clock::now();
        sys_end = system_clock::now();
        tt = system_clock::to_time_t ( sys_begin );
        std::cout << "Program began at:\n " << ctime(&tt);
        tt = system_clock::to_time_t ( sys_end );
        std::cout << "Program finished at:\n " << ctime(&tt);
        duration<double> time_span = duration_cast< duration<double> >(end - begin);
        std::cout << "Running time: " << time_span.count() << " seconds.\n" << std::endl;
        //int x = time_span.count();
        //std::cout << x << std::endl;
    };
};
#endif

