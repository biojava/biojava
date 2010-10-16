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
public class SimpleTranslationTable 
      extends AbstractTranslationTable
      implements Serializable {
  private final Map transMap;
  private final FiniteAlphabet source;
  private final Alphabet target;

  public Alphabet getSourceAlphabet() {
    return source;
  }

  public Alphabet getTargetAlphabet() {
    return target;
  }

  public Symbol doTranslate(Symbol sym) {
    return (Symbol) transMap.get(sym);
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
  }

  /**
   * Create a new translation table that will translate symbols from source to
   * target.
   * <p>
   * The source alphabet must be finite, as a Map object is used to associate
   * a source Symbol with a target Symbol.
   * The target alphabet need not be finite.
   *
   * @param source  the FiniteAlphabet to translate from
   * @param target  the Alphabet to translate into
   */
  public SimpleTranslationTable(FiniteAlphabet source, Alphabet target) {
    this.source = source;
    this.target = target;
    this.transMap = new HashMap();
  }

  /**
   * Create a new translation table that will translate symbols from source to
   * target.
   * <p>
   * The Map transMap should contain keys in the source alphabet with values in
   * the target alphabet. However, this is not currently checked.
   * <p>
   * The source alphabet must be finite, as a Map object is used to associate
   * a source Symbol with a target Symbol.
   * The target alphabet need not be finite.
   *
   * @param source  the FiniteAlphabet to translate from
   * @param target  the Alphabet to translate into
   */
  public SimpleTranslationTable(
    FiniteAlphabet source, Alphabet target, Map transMap
  ) {
    this.source = source;
    this.target = target;
    this.transMap = transMap;
  }
}
