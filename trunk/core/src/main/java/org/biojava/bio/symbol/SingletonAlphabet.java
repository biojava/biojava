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
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.SingletonList;

/**
 * An alphabet that contains a single atomic symbol.
 *
 * @author Matthew Pocock
 */
public class SingletonAlphabet
extends AbstractAlphabet
implements FiniteAlphabet, Serializable {
  private final AtomicSymbol sym;
  private List alphabets;
  
  public SingletonAlphabet(AtomicSymbol sym) {
    this.sym = sym;
  }

  public List getAlphabets() {
    if(this.alphabets == null) {
      this.alphabets = new SingletonList(this);
    }
    return this.alphabets;
  }
  
  protected boolean containsImpl(AtomicSymbol s) {
    return s == sym;
  }
  
  public String getName() {
    return sym.getName() + "-alphabet";
  }
  
  public SymbolTokenization getTokenization(String name)
  throws NoSuchElementException {
    throw new NoSuchElementException(
      "No parsers associated with " + getName() +
      ": " + name
    );
  }
  
  public Iterator iterator() {
    return Collections.singleton(sym).iterator();
  }
  
  public int size() {
    return 1;
  }
  
  public Annotation getAnnotation() {
    return Annotation.EMPTY_ANNOTATION;
  }
  
  public void addSymbolImpl(AtomicSymbol sym) throws IllegalSymbolException {
    throw new IllegalSymbolException(
      "Can't add symbols to alphabet: " + sym.getName() +
      " in " + getName()
    );
  }
  
  public void removeSymbol(Symbol sym) throws IllegalSymbolException {
    throw new IllegalSymbolException(
      "Can't remove symbols from alphabet: " + sym.getName() +
      " in " + getName()
    );
  }
  
  protected AtomicSymbol getSymbolImpl(List symList)
  throws IllegalSymbolException {
    return (AtomicSymbol) symList.get(0);
  }
}
