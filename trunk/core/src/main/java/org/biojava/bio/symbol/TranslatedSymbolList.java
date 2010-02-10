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

import org.biojava.bio.BioError;

/**
 * Provides a 'translated' view of an underlying SymbolList.
 * <p>
 * This class allows you to translate from one alphabet into another, so
 * for example, you could translate from DNA-triplets into amino-acids. You
 * could also translate from DNA-dinucleotide into the 'twist' structural
 * metric, or any other translation that takes your fancy.
 * <p>
 * The actual mapping from source to view Symbol is encapsulated in a
 * TranslationTable object.
 * <p>
 * The TranslatedSymbolList will be the same length as the source, and each
 * Symbol in the view will correspond to a single Symbol in the source.
 *
 * @author Matthew Pocock
 */
class TranslatedSymbolList
extends AbstractSymbolList implements SymbolList {
  /**
   * The source symbol list to translate.
   */
  private final SymbolList source;

  /**
   * The TranslationTable that will be used to translate source->view symbols
   */
  private final TranslationTable transTable;

    /**
*Obtain the translation table associated with this symbol list
*/


  public TranslationTable getTranslationTable() {
    return transTable;
  }

    /**
*Returns the symbol list associated with this translated symbol list.
*/

  public SymbolList getSource() {
    return source;
  }

  public TranslatedSymbolList(SymbolList source, TranslationTable transTable)
  throws IllegalAlphabetException {
    if(transTable.getSourceAlphabet() != source.getAlphabet()) {
      throw new IllegalAlphabetException(
        "The source alphabet and translation table source alphabets don't match: " +
        source.getAlphabet().getName() + " and " +
        transTable.getSourceAlphabet().getName()
      );
    }

    this.source = source;
    this.transTable = transTable;
  }

  public int length() {
    return source.length();
  }

  public Symbol symbolAt(int indx) {
    try {
      return transTable.translate(source.symbolAt(indx));
    } catch (IllegalSymbolException ire) {
      throw new BioError(
        "I thought that I had checked that the translation table was compatible with " +
        "my source, but apparently something has messed up.", ire
      );
    }
  }

  public Alphabet getAlphabet() {
    return transTable.getTargetAlphabet();
  }
}
