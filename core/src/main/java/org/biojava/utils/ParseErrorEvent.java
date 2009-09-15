/*
 * BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 * http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 * http://www.biojava.org
 */

package org.biojava.utils;

import java.util.EventObject;

/**
 *  Event which signals a bad line when parsing a record.
 *
 * @author     Greg Cox
 */

public class ParseErrorEvent extends EventObject
{
	private final String mMessage;

// Constructors
	/**
	 *  Construct a ParseErrorEvent with no other information.
	 *
	 * @param  theSource  The source of the parse error
	 */
	public ParseErrorEvent(Object theSource)
	{
		this(theSource, null);
	}

	/**
	 *
	 * Construct a ParseErrorEvent with a message.
	 *
	 * @param  theSource  The source of the parse error.
	 * @param  theMessage The message.
	 */

	public ParseErrorEvent(Object theSource, String theMessage)
	{
		super(theSource);
		this.mMessage = theMessage;
	}

// Accessor Methods
	/**
	 * Find the message about this event
	 *
	 * @return The message.
	 */
	public String getMessage()
	{
		return mMessage;
	}

	public String toString()
	{
		StringBuffer representation = new StringBuffer(
			super.toString() +
			"[" +
			this.getMessage() +
			"]");
		return representation.substring(0);
  }
}
