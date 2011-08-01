/* @(#)InputParameters.java	1.0 June 2010
 * 
 *  Copyright (c) 2010 Peter Troshin
 * 
 *  JRONN version: 3.1     
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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import compbio.ronn.ORonn.ResultLayout;

/**
 * Holder for input parameters of the {@link ORonn} class
 * 
 * @author Petr Troshin
 * @version 1.0
 */
final class InputParameters {
    static final String inputKey = "-i=";
    static final String outputKey = "-o=";
    static final String statKey = "-s=";
    static final String disorderKey = "-d=";
    static final String formatKey = "-f=";
    static final String threadKey = "-n=";

    private File input;
    private File output;
    private File statistics;
    private ResultLayout format;
    private float disorder;
    private int threadNum;
    private volatile PrintWriter outWriter;
    private volatile PrintWriter statWriter;

    InputParameters() {
	threadNum = Runtime.getRuntime().availableProcessors();
	disorder = RonnConstraint.DEFAULT_DISORDER;
    }

    File setFilePrm(String filename, final String key) throws IOException {
	File file = null;
	if (filename == null) {
	    throw new IllegalArgumentException("File name is not provided! ");
	}
	filename = filename.trim();
	if (filename.toLowerCase().startsWith(key)) {
	    file = new File(filename.substring(key.length()));
	    if (!file.exists()) {
		file.createNewFile();
	    }
	    if (key.equals(InputParameters.inputKey)) {
		input = file;
	    }
	    if (key.equals(InputParameters.outputKey)) {
		output = file;
	    }
	    if (key.equals(InputParameters.statKey)) {
		statistics = file;
	    }
	}
	return file;
    }

    PrintWriter getOutputWriter() throws FileNotFoundException {
	if (outWriter == null) {
	    synchronized (this) {
		if (outWriter == null) {
		    if (output != null) {
			outWriter = new PrintWriter(output);
		    } else {
			outWriter = new PrintWriter(System.out);
		    }
		}
	    }
	}
	return outWriter;
    }

    ResultLayout parseFormat(String format) {
	// default layout
	if (format == null) {
	    return ResultLayout.VERTICAL;
	}
	format = format.trim().substring(InputParameters.formatKey.length());
	if (format.toUpperCase().equals("V")) {
	    return ResultLayout.VERTICAL;
	}
	if (format.toUpperCase().equals("H")) {
	    return ResultLayout.HORIZONTAL;
	}
	throw new IllegalArgumentException("Unrecognised format: '" + format
		+ "' Output format for results must be either H or V.");
    }

    /**
     * 
     * @param key
     * @return
     * @throws NumberFormatException
     */
    int parseThreadNum(String key) {
	key = key.trim().substring(InputParameters.threadKey.length());
	final int nthreads = Integer.parseInt(key);
	return nthreads;
    }

    @Override
    public String toString() {
	final String def = " [default]";
	final String newLine = System.getProperty("line.separator");
	String value = "input=" + input + newLine;

	if (output != null) {
	    value += "output=" + output + newLine;
	} else {
	    value += "output= standard out" + newLine;
	}
	if (statistics != null) {
	    value += "statistics=" + statistics + newLine;
	} else {
	    value += "statistics= not kept" + newLine;
	}
	value += "disorder=" + getDisorder();
	if (disorder == RonnConstraint.DEFAULT_DISORDER) {
	    value += def;
	}
	value += newLine;
	value += "thread number=" + getThreadNum() + newLine;
	value += "format=" + getFormat();
	if (format == null) {
	    value += def;
	}
	// value += newLine ;
	return value;
    }

    ResultLayout getFormat() {
	return format == null ? ResultLayout.VERTICAL : format;
    }

    float getDisorder() {
	return disorder;
    }

    void setDisorder(final String prm) {
	final float disorder = Float.parseFloat(prm
		.substring(InputParameters.disorderKey.length()));
	// do not set disorder values if it is invalid
	// default disorder set in the constructor
	if ((disorder > 0) || (disorder <= 1)) {
	    this.disorder = disorder;
	}
    }

    void setFormat(final String prm) {
	format = parseFormat(prm);
    }

    void setThreadNum(final String prm) {
	// default number of threads are set in the constructor and
	// equals to the number of cores
	final int procNum = threadNum;
	final int tnum = parseThreadNum(prm);
	if (tnum != 0) {
	    if ((tnum < 1) || (tnum > procNum * 2)) {
		throw new IllegalArgumentException(
			"Number of threads must be between 1 "
				+ "and 2 x number of cores available on the executing machine");
	    }
	    threadNum = tnum;
	}
    }

    private void initStatWriter() throws IOException {
	if (statistics != null) {
	    statWriter = new PrintWriter(new FileWriter(statistics), true);
	} else {
	    statWriter = new PrintWriter(new NullOutputStream());
	}
    }

    public PrintWriter getStatWriter() throws IOException {
	if (statWriter == null) {
	    synchronized (this) {
		if (statWriter == null) {
		    initStatWriter();
		}
	    }
	}
	return statWriter;
    }

    File getInput() {
	return input;
    }

    int getThreadNum() {
	return threadNum;
    }

}