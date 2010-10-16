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
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * A no-frills implementation of a translation table that
 * maps between two alphabets.  The mapping can be either
 * one-to-one or many-to-one.
 *
 * @author David Huen
 */
public class SimpleManyToOneTranslationTable
      extends AbstractManyToOneTranslationTable
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

  protected Set doUntranslate(Symbol sym) {
    return (Set) revMap.get(sym);
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

    // when putting in a association in a many-to-one
    // you will need to verify that an existing relation does
    // not exist: otherwise you will need to remove it prior
    // to putting in the new one.

    // do I have an existing association for this symbol?
    Symbol prevTransSymbol = doTranslate(from);

    if (prevTransSymbol != null) {
        // there is an association, I need to
        // remove the previous target Symbol
        // from the unTranslate mapping
        // before introducing new mapping.

        Set prevUntransSet = (Set) revMap.get(prevTransSymbol);
        prevUntransSet.remove(from);
    }

    // now enter in the new association
    Set sourceSet = (Set) revMap.get(to);

    if (sourceSet == null) {
        // create new set
        sourceSet = new HashSet();
        revMap.put(to, sourceSet);
    }

    sourceSet.add(from);

    transMap.put(from, to);

  }

  /**
   * Construct a new translation table.
   *
   * @param source  the source FiniteAlphabet
   * @param target the target FiniteAlphabet
   */
  public SimpleManyToOneTranslationTable(FiniteAlphabet source, FiniteAlphabet target)
  {

    this.source = source;
    this.target = target;
    this.transMap = new HashMap();
    this.revMap = new HashMap();
  }
}
