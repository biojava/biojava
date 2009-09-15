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
 * Source for a ParseErrorEvent.  A ParseErrorEvent signals a problem parsing
 * a file.
 *
 * @see        ParseErrorEvent
 * @author     Greg Cox
 */

public interface ParseErrorSource
{
    /**
     * Adds a parse error listener to the list of listeners.
     *
     * @param theListener Listener to be added.
     */
    public void addParseErrorListener(ParseErrorListener theListener);

    /**
     * Removes a parse error listener from the list of listeners.
     *
     * @param theListener Listener to be removed.
     */
    public void removeParseErrorListener(ParseErrorListener theListener);
}
