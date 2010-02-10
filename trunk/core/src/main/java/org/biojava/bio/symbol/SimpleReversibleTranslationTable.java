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

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * A no-frills implementation of TranslationTable that uses a Map to map from
 * symbols in a finite source alphabet into a target alphabet.
 *
 * @author Matthew Pocock
 * @author David Huen (refactoring)
 */
public class SimpleReversibleTranslationTable 
      extends AbstractReversibleTranslationTable
      implements Serializable {
  private final Map transMap;
  private final Map revMap;
  private final FiniteAlphabet source;
  private final Alphabet target;

  public Alphabet getSourceAlphabet() {
    return source;
  }

  public Alphabet getTargetAlphabet() {
    return target;
  }

  protected Symbol doTranslate(Symbol sym) {
    return (Symbol) transMap.get(sym);
  }

  protected Symbol doUntranslate(Symbol sym) {
    return (Symbol) revMap.get(sym);
  }

  /**
   * Alter the translation mapping.
   *
   * @param from source AtomicSymbol
   * @param to   target AtomicSymbol to be returned by translate(from)
   * @throws IllegalSymbolException if either from is not in the source
   *         alphabet or to is not in the target alphabet
   */
  public void setTranslation(AtomicSymbol from, AtomicSymbol to)
  throws IllegalSymbolException {
    source.validate(from);
    target.validate(to);
    transMap.put(from, to);
    revMap.put(to, from);
  }

  /**
   * Construct a new translation table.
   *
   * @param source  the source FiniteAlphabet
   * @param target the target FiniteAlphabet
   * @throws IllegalAlphabetException if the alphabets are of different sizes
   */
  public SimpleReversibleTranslationTable(FiniteAlphabet source, FiniteAlphabet target) 
    throws IllegalAlphabetException
  {

    if(source.size() != target.size()) {
      throw new IllegalAlphabetException(
        "Couldn't create translation table as " +
        "the alphabets were different sizes: " +
        source.size() + ":" + source.getName() +
        target.size() + ":" + target.getName()
      );
    }

    this.source = source;
    this.target = target;
    this.transMap = new HashMap();
    this.revMap = new HashMap();
  }
}
