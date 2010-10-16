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
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.biojava.bio.Annotatable;
import org.biojava.bio.Annotation;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.SingletonList;

/**
 * A simple no-frills implementation of the FiniteAlphabet interface.
 *
 * @serial WARNING the serialized version of this class may not be compatible with future versions
 * @author Matthew Pocock
 */
public class SimpleAlphabet
extends AbstractAlphabet
implements Serializable {
  private String name;
  private Annotation annotation;
  private final Set symbols;
  protected transient ChangeForwarder annotationForwarder;

    //BE SURE TO CHANGE THIS VALUE IF YOU CHANGE THE IMPLEMENTATION
    //SUCH THAT SERIALIZATION WILL FAIL.
  private static final long serialVersionUID = 216254146;

  /**
   * A list of alphabets that make up this one - a singleton list containing
   * this.
   */
  private List alphabets;

  public Iterator iterator() {
    return symbols.iterator();
  }

  public String getName() {
    return name;
  }

  /**
   * Assign a name to the alphabet
   * @param name the name you wish to give this alphabet
   */
  public void setName(String name) {
    this.name = name;
  }

  public Annotation getAnnotation() {
    if(annotation == null)
      annotation = new SimpleAnnotation();
    return annotation;
  }

  public int size() {
    return symbols.size();
  }

  protected boolean containsImpl(AtomicSymbol s) {
    return symbols.contains(s);
  }

  protected void addSymbolImpl(AtomicSymbol s)
  throws IllegalSymbolException, ChangeVetoException {
    symbols.add(s);
  }

  public void removeSymbol(Symbol s)
  throws IllegalSymbolException {
    validate(s);
    //change checking should probably happen in AbstractAlphabet but oh well.
    if(hasListeners(Alphabet.SYMBOLS)) {
      ChangeSupport cs = getChangeSupport(Alphabet.SYMBOLS);
      synchronized(cs) {
        ChangeEvent ce = new ChangeEvent(this, Alphabet.SYMBOLS, null, s);
        cs.firePreChangeEvent(ce);
        _removeSymbol(s);
        cs.firePostChangeEvent(ce);
      }
    }else{
        _removeSymbol(s);
    }
  }

    private void _removeSymbol(final Symbol s) {
        
        if(s instanceof AtomicSymbol) {
          symbols.remove(s);
        } else {
          FiniteAlphabet sa = (FiniteAlphabet) s.getMatches();
          Iterator i = sa.iterator();
          while(i.hasNext()) {
            Symbol sym = (Symbol) i.next();
            symbols.remove(sym);
          }
        }
    }

  public List getAlphabets() {
    if(this.alphabets == null) {
      this.alphabets = new SingletonList(this);
    }
    return this.alphabets;
  }

  protected AtomicSymbol getSymbolImpl(List symL)
  throws IllegalSymbolException {
    AtomicSymbol s = (AtomicSymbol) symL.get(0);
    return s;
  }

  public SimpleAlphabet() {
    this(new HashSet(), null);
  }

  public SimpleAlphabet(Set symbols) {
    this(symbols, null);
  }

  public SimpleAlphabet(String name) {
    this(new HashSet(), name);
  }

  public SimpleAlphabet(Set symbols, String name) {
    this.symbols = new HashSet();
    this.name = name;
    this.alphabets = null;

    // this costs, but I am tracking down a ClassCast exception.
    // roll on parameterised types.
    for(Iterator i = symbols.iterator(); i.hasNext(); ) {
      AtomicSymbol a = (AtomicSymbol) i.next();
      this.symbols.add(a);
    }
  }

  protected ChangeSupport getChangeSupport(ChangeType ct){
    ChangeSupport cs = super.getChangeSupport(ct);

    if(annotationForwarder == null &&
      (ct.isMatchingType(Annotatable.ANNOTATION) || Annotatable.ANNOTATION.isMatchingType(ct)))
    {
      annotationForwarder =
              new ChangeForwarder.Retyper(this, cs, Annotation.PROPERTY);
      getAnnotation().addChangeListener(
          annotationForwarder,
          Annotatable.ANNOTATION);
    }
    return cs;
  }

}
