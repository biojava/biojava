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

package org.biojava.bio.program.tagvalue;

/**
 * <code>ParserListener</code> is an immutable pairing of a parser and
 * listener.
 *
 * @author Matthew Pocock
 * @author Keith James
 */
public class ParserListener {
  private final TagValueParser parser;
  private final TagValueListener listener;
  
    /**
     * Creates a new <code>ParserListener</code> instance.
     *
     * @param parser a <code>TagValueParser</code>.
     * @param listener a <code>TagValueListener</code>.
     */
    public ParserListener(TagValueParser parser, TagValueListener listener) {
    this.parser = parser;
    this.listener = listener;
  }
  
    /**
     * <code>getParser</code> returns the parser of the pair.
     *
     * @return a <code>TagValueParser</code>.
     */
    public TagValueParser getParser() {
    return parser;
  }
  
    /**
     * <code>getListener</code> returns the listener of the pair.
     *
     * @return a <code>TagValueListener</code>.
     */
    public TagValueListener getListener() {
    return listener;
  }
}
