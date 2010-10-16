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
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ListTools;

/**
 * Cross product of a list of arbitrary alphabets.  This is the
 * most flexible implementation of CrossProductAlphabet, but it
 * is likely to be possible to produce more efficient implementations
 * for specific tasks.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Greg Cox
 */

class SimpleCrossProductAlphabet
extends AbstractAlphabet
implements Serializable {
  private final Alphabet parent;
  private final List alphas;
  private final HashMap ourSymbols;

  /**
   * Create a cross-product alphabet over the list of alphabets in 'a'.
   */
  public SimpleCrossProductAlphabet(List a)
  throws IllegalAlphabetException {
    this(a, null);
  }

  public SimpleCrossProductAlphabet(List a, Alphabet parent)
  throws IllegalAlphabetException {
    if(a.size() == 0) {
      throw new IllegalAlphabetException(
        "Can't create alphabet for empty list. Use Alphabet.EMPTY_ALPHABET"
      );
    }

    this.parent = parent;
    for(Iterator i = a.iterator(); i.hasNext(); ) {
      Alphabet aa = (Alphabet) i.next();
      if(! (aa instanceof FiniteAlphabet) ) {
        throw new IllegalAlphabetException(
          "Can't create a SimpleCrossProductAlphabet over non-fininte alphabet " +
          aa.getName() + " of type " + aa.getClass()
        );
      }
    }
    alphas = ListTools.createList(a);
    ourSymbols = new HashMap();
    populateSymbols(new ArrayList());
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

  public Iterator iterator() {
    return ourSymbols.values().iterator();
  }

  private void populateSymbols(List symList) {
    if (symList.size() == alphas.size()) {
      putSymbol(symList);
    } else {
      int indx = symList.size();
      FiniteAlphabet a = (FiniteAlphabet) alphas.get(indx);
      Iterator i = a.iterator();
      if(i.hasNext()) {
        symList.add(i.next());
        populateSymbols(symList);
        while (i.hasNext()) {
          symList.set(indx, i.next());
          populateSymbols(symList);
        }
        symList.remove(indx);
      }
    }
  }

  private AtomicSymbol putSymbol(List s) {
    AtomicSymbol ss;
    if(parent != null) {
      try {
        ss = (AtomicSymbol) parent.getSymbol(s);
      } catch (IllegalSymbolException ise) {
        throw new BioError("Balls up - couldn't fetch symbol from parent", ise);
      }
    } else {
      try {
        ss = (AtomicSymbol) AlphabetManager.createSymbol(
          Annotation.EMPTY_ANNOTATION, s, this
        );
      } catch (IllegalSymbolException ise) {
        throw new BioError(

          "Assertion Failure: Should have a legal symbol: " + s, ise
        );
      }
    }
    ourSymbols.put(ss.getSymbols(), ss);
    return ss;
  }

  protected boolean containsImpl(AtomicSymbol s) {
    return ourSymbols.values().contains(s);
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
    return ourSymbols.size();
  }

  public Annotation getAnnotation() {
    return Annotation.EMPTY_ANNOTATION;
  }

  public List getAlphabets() {
    return alphas;
  }

  protected AtomicSymbol getSymbolImpl(List l)
  throws IllegalSymbolException {
    AtomicSymbol cps;
    cps = (AtomicSymbol) ourSymbols.get(l);

    if(cps == null) {
      throw new IllegalSymbolException(
        "Can't find symbol for " + l +
        " in alphabet " + getName()
      );
    }

    return cps;
  }

  protected void addSymbolImpl(AtomicSymbol sym)
  throws IllegalSymbolException {
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
}
