// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
// Created by Francois Nedelec on 27/9/2008.

#include "time_date.h"

#include <cstring>
#include <sys/time.h>

//------------------------------------------------------------------------------
#pragma mark - Wall time

/**
 This get current time from the C-library functions time() and ctime_r()
 */
void TimeDate::get_date(char * buf, size_t buf_size)
{
    if ( buf_size > 25 )
    {
        time_t now = time(nullptr);
        //asctime_r(localtime(&now), buf);
#if ( 0 )
        ctime_r(&now, buf);
#else
        strncpy(buf, ctime(&now), buf_size);
#endif
        // remove new line:
        buf[24] = 0;
    }
    else if ( buf_size > 0 )
        buf[0] = 0;
    // terminate string:
    if ( buf_size > 0 )
        buf[buf_size-1] = 0;
}


void TimeDate::get_date(char * buf, size_t buf_size, bool no_year)
{
    if ( buf_size > 25 )
    {
        TimeDate::get_date(buf, buf_size);
        // remove year:
        if ( no_year )
            buf[19] = 0;
    }
}


char const* TimeDate::date_string()
{
    static char buf[32];
    get_date(buf, sizeof(buf));
    return buf;
}


int TimeDate::days_since_2000()
{
    time_t now = time(nullptr);
    tm * loc = localtime(&now);
    return loc->tm_year - 2000 + loc->tm_yday;
}


/** Attention: this will fail after Jan. 2038 */
time_t TimeDate::seconds_since_1970()
{
    return time(nullptr);
}

/** Attention: this will fail after Jan. 2038 */
time_t TimeDate::seconds_since_2000()
{
    return time(nullptr) - 946684800;
}


int TimeDate::year()
{
    time_t now = time(nullptr);
    tm * loc = localtime(&now);
    return loc->tm_year;
}


int TimeDate::day_of_the_year()
{
    time_t now = time(nullptr);
    tm * loc = localtime(&now);
    return loc->tm_yday;
}


int TimeDate::hours_today()
{
    time_t now = time(nullptr);
    tm * loc = localtime(&now);
    return loc->tm_hour;
}


double TimeDate::seconds_today()
{
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return (double)tv.tv_sec + 1e-6 * tv.tv_usec;
}


double TimeDate::milliseconds()
{
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return 1000 * tv.tv_sec + tv.tv_usec / 1000.0;
}

double TimeDate::microseconds()
{
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return tv.tv_usec;
}


//------------------------------------------------------------------------------
#pragma mark - CPU time, from the C-library function clock()


double TimeDate::processor_time(clock_t& clk)
{
    clock_t now = clock();
    double sec = double( now - clk ) / CLOCKS_PER_SEC;
    clk = now;
    return sec;
}

