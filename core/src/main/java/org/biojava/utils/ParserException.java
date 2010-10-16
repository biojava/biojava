/*
 *                    BioJava development code
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

package org.biojava.utils;

/**
 * Exception thrown when an error occurs in document parsing.
 * It may optionally include the following fields:
 *
 * <pre>
 *   Locator (file name, URL, etc.)
 *   Line number (negative for unknown)
 *   The text of the actual offending line (Null if unknown)
 *   Character offset (negative for unknown)
 * </pre>
 * @author Matthew Pocock
 * @author Greg Cox
 */

public class ParserException extends Exception {
    private String locator = null;
    private int lineNumber = -1;
    private int character = -1;
    private String line = null;

    public ParserException(String detail) {
	super(detail);
    }

    public ParserException(String detail, String locator) {
	super(detail);
	this.locator = locator;
    }

    public ParserException(String detail, String locator, int line) {
	super(detail);
	this.locator = locator;
	this.lineNumber = line;
    }

    public ParserException(String detail,
			  String locator,
			  int lineNumber,
			  String line)
    {
	super(detail);
	this.locator = locator;
	this.lineNumber = lineNumber;
	this.line = line;
    }

    public ParserException(String detail,
			  String locator,
			  int lineNumber,
			  String line,
			  int character)
    {
	super(detail);
	this.locator = locator;
	this.lineNumber = lineNumber;
	this.character = character;
	this.line = line;
    }

    public ParserException(Throwable t) {
      super(t);
    }

    /**
     * @deprecated use new ParserException(detail, t)
     */
    public ParserException(Throwable t, String detail) {
      this(detail, t);
    }
    
    public ParserException(String message, Throwable cause) {
      super(message, cause);
    }

    public ParserException(Throwable t, String detail, String locator) {
	super(detail, t);
	this.locator = locator;
    }

    public ParserException(Throwable t, String detail, String locator, int line) {
	super(detail, t);
	this.locator = locator;
	this.lineNumber = line;
    }

    public ParserException(
        Throwable t,
        String detail,
			  String locator,
			  int lineNumber,
			  String line)
    {
	super(detail, t);
	this.locator = locator;
	this.lineNumber = lineNumber;
	this.line = line;
    }

    public ParserException(
        Throwable t,
        String detail,
			  String locator,
			  int lineNumber,
			  String line,
			  int character)
    {
	super(detail, t);
	this.locator = locator;
	this.lineNumber = lineNumber;
	this.character = character;
	this.line = line;
    }

    /**
     * Get a locator for the stream which caused this exception.
     *
     * @return A locator string, or <code>null</code> if none is
     *         known.
     */

    public String getLocator() {
	return locator;
    }

    /**
     * Get the line number in the stream where this exception occured.
     *
     * @return A positive integer line number, or -1 if not known.
     */

    public int getLineNumber() {
	return lineNumber;
    }

    /**
     * Get the character offset in the line where an error was detected.
     *
     * @return The first character in the line where the parser detected
     *         an error, or -1 if the exception effects the whole line.
     */

    public int getCharacterOffset() {
	return character;
    }

    /**
     * Get the text of the line where the exception occured.
     *
     * @return The text of the line, or <code>null</code> if not known.
     */

    public String getLine() {
	return line;
    }

    /**
     * Represent this exception as a string.  This includes
     * the default exception toString representation, followed
     * by details of the location where the error occured, if
     * they were supplied when constructing this exception.
     *
     * @return A string representation of this exception.
     */

    public String toString() {
	StringBuffer sb = new StringBuffer(super.toString());
	if (locator != null) {
	    sb.append('\n');
	    sb.append("Parsing location: ");
	    sb.append(locator);
	}

	if (lineNumber >= 0) {
	    sb.append('\n');
	    sb.append("Parsing line: ");
	    sb.append(lineNumber);
	}

	if (line != null) {
	    sb.append('\n');
	    sb.append(line);
	    if (character >= 0) {
		sb.append('\n');
		for (int i = 0; i < character; ++i)
		    sb.append(' ');
		sb.append('^');
	    }
	}
	sb.append('\n');
	return sb.substring(0);
    }
}
