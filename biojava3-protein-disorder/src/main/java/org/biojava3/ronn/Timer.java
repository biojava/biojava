/* 
 * @(#)Timer.java	1.0 June 2010
 * 
 * Copyright (c) 2010 Peter Troshin
 *  
 *        BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava3.ronn;

import java.util.concurrent.TimeUnit;

/**
 * A simple timer, calculates the time interval between two events. Keeps two
 * counters, one for long time intervals, to measure time between the start and
 * end of the application for instance, and another for short events, to measure
 * how long it took to reach a next block of code.
 * 
 * @author Peter Troshin
 * @version 1.0
 * @since 3.0.2
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
