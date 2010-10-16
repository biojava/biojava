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

import java.io.ObjectStreamException;
import java.io.Serializable;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.CrossProductTokenization;
import org.biojava.bio.seq.io.NameTokenization;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;

/**
 * <p>
 * An abstract implementation of <code>Alphabet</code>.
 * </p>
 *
 * <p>
 * This provides the frame-work for maintaining the SymbolParser <-> name
 * mappings and also for the ChangeListeners.
 * </p>
 *
 * <p>
 * This class is for developers to derive from, not for use directly.
 * </p>
 *
 * @author Matthew Pocock
 * @author Greg Cox
 * @author Francois Pepin
 * @author Mark Schreiber
 * @since 1.1
 */
public abstract class AbstractAlphabet
  extends
    AbstractChangeable
  implements
    FiniteAlphabet,
    Serializable
{
  private final Map tokenizationsByName;
  private final Map ambCache;
  public static final long serialVersionUID = -3043128839927615753l;
  
  /*
   * this field records if the alpha has been registered with the AlphabetManager
   * in the current VM. It is important for serialization, if the alpha has been registered
   * it needs to be registered at the other end.
   */
  private boolean registered;
  
  /**
   * Used by the AlphabetManager to inform an Alphabet that is has been
   * registered and that is should be registered if it is transported to a new
   * AlphabetManager on another VM
   */
  void setRegistered(boolean value){
      this.registered = value;
  }
  boolean getRegsitered(){
      return(registered);
  }
  
  
  
  {
    tokenizationsByName = new HashMap();
    ambCache = new HashMap();
    registered = false;
  }


  /**
   * To prevent duplication of a what should be a
   * single instance of an existing alphabet. This method
   * was written as protected so that subclasses even from
   * other packages will inherit it. It should only be overridden
   * with care.
   */
  protected Object readResolve() throws ObjectStreamException{
      
      if(AlphabetManager.registered(this.getName())){
          return AlphabetManager.alphabetForName(this.getName());
      }else{
          if(this.registered){
              AlphabetManager.registerAlphabet(this.getName(), this);
          }
          return this;
      }
  }

  /**
   * <p>
   * Assigns a symbol parser to a String object.
   * </p>
   *
   * <p>
   * Afterwards, the parser can be retrieved using the
   * getTokenization(Sting name) method.
   * </p>
   */
  public void putTokenization(String name, SymbolTokenization parser) {
    tokenizationsByName.put(name, parser);
  }

  public SymbolTokenization getTokenization(String name)
        throws NoSuchElementException, BioException
    {
        SymbolTokenization toke = (SymbolTokenization) tokenizationsByName.get(name);
        if(toke == null) {
            if(name.equals("name")) {
                if (getAlphabets().size() == 1) {
                    toke = new NameTokenization(this);
                } else {
                    toke = new CrossProductTokenization(this);
                }
                putTokenization(name, toke);
            } else if (name.equals("default")) {

              if (tokenizationsByName.containsKey("token"))
                toke= (SymbolTokenization)tokenizationsByName.get("token");
              else
                toke= (SymbolTokenization)getTokenization("name");
              putTokenization(name, toke);

            }
            else
            {
                throw new NoSuchElementException("There is no tokenization '" + name +
                                                 "' defined in alphabet " + getName());
            }
        }
        return toke;
    }

    public final Symbol getAmbiguity(Set syms)
        throws IllegalSymbolException
    {
        if (syms.size() == 0) {
            return getGapSymbol();
        } else if (syms.size() == 1) {
            Symbol sym = (Symbol) syms.iterator().next();
            validate(sym);
            return sym;
        } else {
            Symbol s = (Symbol) ambCache.get(syms);
            if(s == null) {
                for (Iterator i = syms.iterator(); i.hasNext(); ) {
                    validate((Symbol) i.next());
                }


                s = getAmbiguityImpl(syms);
                ambCache.put(new HashSet(syms), s);
            }
            return s;
        }
    }

    /**
     * Backend for getAmbiguity, called when it is actually necessarly to create a new symbol.
     * By default, calls AlphabetManager.createSymbol.
     *
     * @since 1.3
     */

    protected Symbol getAmbiguityImpl(Set syms)
        throws IllegalSymbolException
    {
        return AlphabetManager.createSymbol(
                    Annotation.EMPTY_ANNOTATION,
                    syms, this
                );
    }

  public final Symbol getSymbol(List syms)
  throws IllegalSymbolException {
    if (syms.size() == 1) {
      Symbol s = (Symbol) syms.get(0);
      validate(s);
      return s;
    }

    List alphas = getAlphabets();

    if(alphas.size() != syms.size()) {
      throw new IllegalSymbolException(
        "Can't retrieve symbol as symbol list is the wrong length " +
        syms.size() + ":" + alphas.size()
      );
    }

    Iterator si = syms.iterator();
    int atomic = 0;
    while(si.hasNext()) {
      Symbol s = (Symbol) si.next();
      //Alphabet a = (Alphabet) ai.next();
      //a.validate(s); // very expensive for requent fetches!
      if(s instanceof AtomicSymbol) {
        atomic++;
      }
    }

    if(atomic == syms.size()) {
      return getSymbolImpl(syms);
    } else {
      return AlphabetManager.createSymbol(
        Annotation.EMPTY_ANNOTATION,
        syms, this
      );
    }
  }

  protected abstract AtomicSymbol getSymbolImpl(List symList)
  throws IllegalSymbolException;

  protected abstract void addSymbolImpl(AtomicSymbol s)
  throws IllegalSymbolException, ChangeVetoException;

  public final void addSymbol(Symbol s)
  throws IllegalSymbolException, ChangeVetoException {
    if(s == null) {
      throw new IllegalSymbolException(
        "You can not add null as a symbol to alphabet " + getName()
      );
    }

    if(hasListeners()) {
      ChangeSupport cs = getChangeSupport(Alphabet.SYMBOLS);
      synchronized(cs) {
        ChangeEvent ce = new ChangeEvent(this, Alphabet.SYMBOLS, s, null);
        cs.firePreChangeEvent(ce);
        doAddSymbol(s);
        cs.firePostChangeEvent(ce);
      }
    } else {
      doAddSymbol(s);
    }
  }

  private void doAddSymbol(Symbol s)
  throws IllegalSymbolException, ChangeVetoException {
    Alphabet sa = s.getMatches();
    if(!(sa instanceof FiniteAlphabet)) {
      throw new IllegalSymbolException(
        "Can't add symbol " + s.getName() +
        " as it matches an infinite number of symbols."
      );
    } else {
      for(Iterator si = ((FiniteAlphabet) sa).iterator(); si.hasNext(); ) {
        addSymbolImpl((AtomicSymbol) si.next());
      }
    }
  }

  public final boolean contains(Symbol sym) {
    if(sym instanceof AtomicSymbol) {
      return containsImpl((AtomicSymbol) sym);
    } else {
      if(sym == null) {
        throw new NullPointerException("Symbols can't be null");
      }
      FiniteAlphabet matches = (FiniteAlphabet) sym.getMatches();
      if(matches.size() == 0) {
        //System.out.println("Got empty symbol " + sym.getName());
        if(sym.equals(AlphabetManager.getGapSymbol())) {
          //System.out.println("Global gap symbol");
          return true;
        } else if(sym instanceof BasisSymbol) {
          if(((BasisSymbol) sym).getSymbols().size() == getAlphabets().size()) {
            //System.out.println("Basis symbol and the right length");
            return true;
          }
        }
        //System.out.println("Empty symbol and not basis - let's accept it.");
        return true;
      }
      for(Iterator i = matches.iterator(); i.hasNext(); ) {
        AtomicSymbol s = (AtomicSymbol) i.next();
        if(!containsImpl(s)) {
          return false;
        }
      }
      return true;
    }
  }

  public final Symbol getGapSymbol() {
    return AlphabetManager.getGapSymbol(getAlphabets());
  }

  public final void validate(Symbol sym)
  throws IllegalSymbolException {
    if(!contains(sym)) {
      throw new IllegalSymbolException(
        "Symbol " + sym.getName() + " not found in alphabet " + this.getName()
      );
    }
  }

  protected abstract boolean containsImpl(AtomicSymbol s);

  /*

  public boolean equals(Object o) {
    if(o == this) {
      return true;
    }

    if(!(o instanceof FiniteAlphabet)) {
      return false;
    }

    FiniteAlphabet that = (FiniteAlphabet) o;

    if(this.size() != that.size()) {
      return false;
    }

    for(Iterator i = that.iterator(); i.hasNext(); ) {
      if(!this.contains((AtomicSymbol) i.next())) {
        return false;
      }
    }

    return true;
  }

  */

  public String toString() {
    return getName();
  }

  protected AbstractAlphabet() {}
}

