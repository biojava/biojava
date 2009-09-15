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

package org.biojava.bio.symbol;

import java.util.Set;

/**
 * A translation table that will handle the many-to-one mappings
 * that you see, for example, with genetic codes.
 * <p>
 * It differs from a ReversibleTranslationTable in that the reverse
 * translation returns a Set of Symbols in the source alphabet that
 * translate to give that Symbol in the target Alphabet.
 *
 * @author David Huen
 */
public interface ManyToOneTranslationTable extends TranslationTable {
  /**
   * Translate a single symbol from target alphabet to the source alphabet.
   *
   * @param sym the Symbol to reverse-translate (member of target alphabet)
   * @return a Set containing symbols that translate to the specified Symbol.
   * @throws IllegalSymbolException if sym is not a member of the target
   *         alphabet
   */
  public Set untranslate(Symbol sym) throws IllegalSymbolException;
}
