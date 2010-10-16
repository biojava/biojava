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

/**
 * Listener for a ParseErrorEvent. A ParseErrorEvent signals a problem
 * parsing a file.
 *
 * @see        ParseErrorEvent
 * @author     Greg Cox
 */

public interface ParseErrorListener extends java.util.EventListener
{
    /**
     * Method called when the parser encounters a bad line.
     *
     * @param theEvent The event that contains the line and token.
     */
    public abstract void BadLineParsed(ParseErrorEvent theEvent);
}
