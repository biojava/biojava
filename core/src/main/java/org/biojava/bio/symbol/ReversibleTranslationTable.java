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

/**
 * A translation table that can also translate from the target to source
 * alphabet.
 * <p>
 * I guess this is encapsulates an invertible function, and the untranslate
 * method is the inverse operation to translate.
 * <p>
 * It is assumed that untranslate(translate(x)) = x for all x in the source
 * alphabet, and that translate(untranslate(y)) = y for all y in the target
 * alphabet. Note, one interesting sub-set of reversible transforms are of the
 * form translate(x) = untranslate(x), and represent 'mirror image'
 * transformations.
 *
 * @author Matthew Pocock
 */
public interface ReversibleTranslationTable extends TranslationTable {
  /**
   * Translate a single symbol from target alphabet to the source alphabet.
   *
   * @param sym the Symbol to translate (member of target alphabet)
   * @return the translated version of sym (member of source alphabet)
   * @throws IllegalSymbolException if sym is not a member of the target
   *         alphabet
   */
  public Symbol untranslate(Symbol sym) throws IllegalSymbolException;
}
