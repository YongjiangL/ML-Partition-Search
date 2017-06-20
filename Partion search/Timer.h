// ***************************************************************
//  timerec   version:  1.0   ¡¤  date: 01 - 20 - 2007
//  -------------------------------------------------------------
//  
//  -------------------------------------------------------------
//  Copyright (C) 2007 - Tsinghua University. All Rights Reserved
// ***************************************************************

/********************************************************************
	created:	2007/01/20
	created:	20:1:2007   0:21
	filename: 	e:\projects\C++\ContrastGraph\src\util\timerec.h
	file path:	e:\projects\C++\ContrastGraph\src\util
	file base:	timerec
	file ext:	h
	author:		clipse
	
	purpose:	Time Recorder Class
*********************************************************************/

#ifndef _TIMEREC_H_
#define _TIMEREC_H_ 1

#ifdef _WIN32
#pragma once
#pragma warning(disable:4996)
#pragma warning(disable:4018)
#pragma warning(disable:4267)
#endif

#include <vector>
using namespace std;

#ifdef _WIN32
#include <sys/TIMEB.h>
#else
#include <sys/time.h>
#endif

#include <assert.h>

namespace TSINGHUA_CLIPSE_UTIL{

	class TimeRecorder
	{
	public:
		  TimeRecorder ();	// set the current time as the starttime
		  ~TimeRecorder ();

	#ifdef _WIN32
		  TimeRecorder (_timeb *initTime);	// set initTime as the starttime
	#else
		  TimeRecorder (struct timeval *initTime);	// set initTime as the starttime
	#endif

		void check ();		// set a check point
		int getCheckNum ();	// get the number of check points
		double diffTime (int index1, int index2);	// return the time span between check point 'index1' and 'index2', base on 1
		double diffTime (int span);	// return the timespan between the checkpoint 'span'  and starttime

	private:

	#ifdef _WIN32
		_timeb st;
		vector <_timeb>checkPoints;
		double diff_timeval (_timeb *a, _timeb *b);	//return the timespan as double
	#else
		struct timeval st;
		vector < struct timeval >checkPoints;
		double diff_timeval (struct timeval *a, struct timeval *b);	//return the timespan as double
	#endif
	};

}

#endif
