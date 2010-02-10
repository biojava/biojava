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

/**
 * @author Primary author unknown
 * @author Greg Cox
 */

package org.biojava.bio.symbol;

import java.io.Serializable;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.CrossProductTokenization;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.Unchangeable;

class InfiniteCrossProductAlphabet
  extends
    Unchangeable
  implements
    Alphabet,
    Serializable
{
  private final List alphas;

  InfiniteCrossProductAlphabet(List alphas) {
    this.alphas = alphas;
  }

  public Symbol getAmbiguity(Set symSet) {
    throw new BioError("Not implemented yet");
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

  public boolean contains(Symbol s) {
    if(! (s instanceof AtomicSymbol)) {
      Alphabet ma = s.getMatches();
      if(ma instanceof FiniteAlphabet) {
        for(Iterator i = ((FiniteAlphabet) ma).iterator(); i.hasNext(); ) {
          if(!contains((Symbol) i.next())) {
            return false;
          }
        }
        return true;
      } else {
        throw new BioError(
          "Problem: Can't work out if I contain ambiguity symbol " + s.getName()
        );
      }
    } else {
      AtomicSymbol cs = (AtomicSymbol) s;

      List sl = cs.getSymbols();
      if(sl.size() != alphas.size()) {
        return false;
      }

      Iterator ai = alphas.iterator();
      Iterator si = sl.iterator();

      while(ai.hasNext() && si.hasNext()) {
        Alphabet aa = (Alphabet) ai.next();
        Symbol ss = (Symbol) si.next();
        if(!aa.contains(ss)) {
          return false;
        }
      }

      return true;
    }
  }

  public void validate(Symbol s) throws IllegalSymbolException {
    if(!this.contains(s)) {
	    throw new IllegalSymbolException(
        "CrossProductAlphabet " + getName() + " does not accept " + s.getName() +
        " as it is not an instance of CrossProductSymbol or " +
        " an AmbiguitySymbol over a subset of symbols in this alphabet."
      );
    }
  }

  public Annotation getAnnotation() {
    return Annotation.EMPTY_ANNOTATION;
  }

  public List getAlphabets() {
    return alphas;
  }

  public Symbol getSymbol(List sList)
  throws IllegalSymbolException {
    if(sList.size() != alphas.size()) {
      throw new IllegalSymbolException(
        "List of symbols is the wrong length (" + alphas.size() +
        ":" + sList.size() + ")"
      );
    }

    Iterator ai = alphas.iterator();
    Iterator si = sList.iterator();

    while(ai.hasNext() && si.hasNext()) {
      Alphabet aa = (Alphabet) ai.next();
      Symbol ss = (Symbol) si.next();
      if(!aa.contains(ss)) {
        throw new IllegalSymbolException(
          "CrossProductAlphabet " + getName() + " does not accept " + sList +
          " as symbol " + ss.getName() + " is not a member of the alphabet " +
          aa.getName()
        );
      }
    }

    return AlphabetManager.createSymbol(Annotation.EMPTY_ANNOTATION, sList, this);
  }

  public Symbol getGapSymbol() {
    return AlphabetManager.getGapSymbol(getAlphabets());
  }

  public SymbolTokenization getTokenization(String name)
  throws NoSuchElementException, BioException {
      if(name == "name") {
          return new CrossProductTokenization(this);
      }
    throw new NoSuchElementException(
      "No parser for " + name + " is defined for " + getName()
    );
  }
}
