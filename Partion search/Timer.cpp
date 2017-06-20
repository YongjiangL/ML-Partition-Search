#include "stdafx.h"
#include "Timer.h"

#include <assert.h>

#define _MACROSEC_ 1000000

using namespace TSINGHUA_CLIPSE_UTIL;

TimeRecorder::TimeRecorder ()
{
#ifdef _WIN32
	_ftime(&st);
#else
	gettimeofday (&st, 0);
#endif
	checkPoints.push_back (st);
}

#ifdef _WIN32
TimeRecorder::TimeRecorder(_timeb *initTime)
{
	checkPoints.push_back (*initTime);
}
#else
TimeRecorder::TimeRecorder (struct timeval *initTime)
{
	checkPoints.push_back (*initTime);
}
#endif

TimeRecorder::~TimeRecorder ()
{
}


#ifdef _WIN32
double
TimeRecorder::diff_timeval (_timeb *a, _timeb *b)
{
	return (double)(b->time-a->time)+(double)(b->millitm-a->millitm)/1000;
}
#else
double
TimeRecorder::diff_timeval (struct timeval *a, struct timeval *b)
{
	assert (a != 0x00 && b != 0x00);
	double la = (double) a->tv_sec + ((double) a->tv_usec) / _MACROSEC_;
	double lb = (double) b->tv_sec + ((double) b->tv_usec) / _MACROSEC_;
	return lb - la;
}
#endif

void
TimeRecorder::check ()
{
#ifdef _WIN32
	_ftime(&st);
#else
	gettimeofday (&st, 0);
#endif
	checkPoints.push_back (st);
}


int
TimeRecorder::getCheckNum ()
{
	return checkPoints.size () - 1;
}

double
TimeRecorder::diffTime (int index1, int index2)
{
	assert (index1 < checkPoints.size () && index2 < checkPoints.size ());
	return diff_timeval (&checkPoints[index1], &checkPoints[index2]);
}

double
TimeRecorder::diffTime (int span)
{
	assert (span < checkPoints.size ());
	return diff_timeval (&checkPoints[0], &checkPoints[span]);
}
