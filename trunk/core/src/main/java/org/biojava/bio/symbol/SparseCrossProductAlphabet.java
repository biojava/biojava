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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.cache.WeakValueHashMap;

/**
 * Cross product of a list of arbitrary alphabets.  This is a memory efficient
 * implementation of CrossProductAlphabet that instantiates symbols as they are
 * needed. This is required as alphabets can get prohibatively large very
 * quickly (e.g. align 200 proteins & you need 20^200 Symbols).
 *
 * @author Matthew Pocock
 * @author Greg Cox
 */

class SparseCrossProductAlphabet
extends AbstractAlphabet
implements Serializable {
  private final int size;
  private final List alphas;
  private final Map knownSymbols;

  SparseCrossProductAlphabet(List alphas) {
    this.alphas = alphas;
    knownSymbols = new WeakValueHashMap();
    int size = 1;
    for(Iterator i = alphas.iterator(); i.hasNext(); ) {
      FiniteAlphabet a = (FiniteAlphabet) i.next();
      size *= a.size();
    }
    this.size = size;
  }

  protected ChangeSupport generateChangeSupport() {
      for (Iterator i = alphas.iterator(); i.hasNext(); ) {
          Alphabet a = (Alphabet) i.next();
          if (!a.isUnchanging(Alphabet.SYMBOLS)) {
              return new ChangeSupport();
          }
      }
      return new ChangeSupport(Collections.singleton(Alphabet.SYMBOLS));
  }

  public String getName() {
    StringBuffer name = new StringBuffer("(");
    for (int i = 0; i < alphas.size(); ++i) {
            Alphabet a = (Alphabet) alphas.get(i);
            name.append(a.getName());
            if (i < alphas.size() - 1) {
        name.append(" x ");
      }
    }
    name.append(")");
    return name.substring(0);
  }

  public int size() {
    return size;
  }

  protected boolean containsImpl(AtomicSymbol s) {
    return knownSymbols.values().contains(s);
  }

  public Annotation getAnnotation() {
    return Annotation.EMPTY_ANNOTATION;
  }

  public List getAlphabets() {
    return alphas;
  }

  public Iterator iterator() {
    return new SparseIterator(this);
  }

  protected AtomicSymbol getSymbolImpl(List sList)
  throws IllegalSymbolException {
    AtomicSymbol s;
    s = (AtomicSymbol) knownSymbols.get(sList);

    if(s == null) {
      Iterator si = sList.iterator();
      Iterator ai = getAlphabets().iterator();
      while(ai.hasNext()) {
        ((Alphabet) ai.next()).validate((Symbol) si.next());
      }

      List l = new ArrayList(sList);
      s = (AtomicSymbol) AlphabetManager.createSymbol(
        Annotation.EMPTY_ANNOTATION, l, this
      );
      knownSymbols.put(s.getSymbols(), s);
    }

    return s;
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

  private static class SparseIterator implements Iterator {
    private Alphabet parent;
    private FiniteAlphabet []alphas;
    private Iterator []symI;
    private AtomicSymbol []as;
    private boolean hasNext;
    private List symList;

    public SparseIterator(FiniteAlphabet parent) {
      this.parent = parent;
      this.alphas = (FiniteAlphabet []) parent.getAlphabets().toArray(new FiniteAlphabet[0]);
      this.symI = new Iterator[this.alphas.length];
      this.as = new AtomicSymbol[this.alphas.length];
      this.hasNext = true;

      for(int i = 0; i < this.alphas.length; i++) {
        this.symI[i] = alphas[i].iterator();
        if(!symI[i].hasNext()) {
          this.hasNext = false;
          return;
        }
        this.as[i] = (AtomicSymbol) symI[i].next();
      }

      symList = Arrays.asList(as);
    }

    public boolean hasNext() {
      return hasNext;
    }

    public Object next() {
      try {
        Symbol sym = parent.getSymbol(symList);

        for(int i = 0; i <= alphas.length; i++) {
          if(i == alphas.length) {
            hasNext = false;
          } else {
            if(!symI[i].hasNext()) {
              symI[i] = alphas[i].iterator();
              as[i] = (AtomicSymbol) symI[i].next();
            } else {
              as[i] = (AtomicSymbol) symI[i].next();
              break;
            }
          }
        }
        return sym;
      } catch (IllegalSymbolException ise) {
        throw new BioError( "Assertion Failure: I should contain this symbol", ise);
      }
    }

    public void remove()
    throws UnsupportedOperationException {
      throw new UnsupportedOperationException();
    }
  }
}
