/* 
 * @(#)Timer.java	1.0 June 2010
 * 
 * Copyright (c) 2010 Peter Troshin
 *  
 * JRONN version: 3.1     
 *   
 *  This library is free software; you can redistribute it and/or modify it under the terms of the
 *  Apache License version 2 as published by the Apache Software Foundation
 * 
 *  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 *  even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache 
 *  License for more details.
 * 
 *  A copy of the license is in apache_license.txt. It is also available here:
 * see: http://www.apache.org/licenses/LICENSE-2.0.txt
 * 
 * Any republication or derived work distributed in source code form
 * must include this copyright and license notice.
 */
package compbio.ronn;

import java.util.concurrent.TimeUnit;

/**
 * A simple timer, calculates the time interval between two events. Keeps two
 * counters, one for long time intervals, to measure time between the start and
 * end of the application for instance, and another for short events, to measure
 * how long it took to reach a next block of code.
 * 
 * @author Petr Troshin
 * @version 1.0
 * 
 */
public class Timer {

    private long checkPoint;
    private final long startTime;
    private TimeUnit reportTimeUnit;

    public Timer() {
	startTime = System.nanoTime();
	checkPoint = startTime;
	// set default time unit for reporting
	reportTimeUnit = TimeUnit.SECONDS;
    }

    public Timer(final TimeUnit reportIn) {
	this();
	reportTimeUnit = reportIn;
    }

    public void checkPoint() {
	checkPoint = System.nanoTime();
    }

    long getStepTime(final TimeUnit tunit) {
	final long duration = tunit.convert(System.nanoTime() - checkPoint,
		TimeUnit.NANOSECONDS);
	checkPoint();
	return duration;
    }

    long getStepTime() {
	return getStepTime(reportTimeUnit);
    }

    long getTotalTime() {
	return getTotalTime(reportTimeUnit);
    }

    long getTotalTime(final TimeUnit tunit) {
	return tunit.convert(System.nanoTime() - startTime,
		TimeUnit.NANOSECONDS);
    }
}
