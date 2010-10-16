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
 * Communication interface between Parser and a TagValueListener that allows
 * listeners to request that a parser/listener pair be pushed onto the stack to
 * handle the current tag.
 *
 * @author Mathew Pocock
 * @since 1.2
 */
public interface TagValueContext {
  /**
   * <p>
   * Push a parser and listener pair onto the parser stack.
   * </p>
   *
   * <p>
   * This will result in the parser using subParser to process all values of the
   * current tag.
   * </p>
   */
  void pushParser(TagValueParser subParser, TagValueListener subListener);
}
