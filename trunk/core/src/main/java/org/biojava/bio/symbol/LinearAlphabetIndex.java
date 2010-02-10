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

import java.lang.ref.Reference;
import java.lang.ref.WeakReference;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeAdapter;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/*
 * Implementatoin of AlphabetIndex that stores the symbols in an array and does
 * a linear shearch through the list for a given symbol to find its index.
 *
 * @author Matthew Pocock
 * @since 1.1
 */
class LinearAlphabetIndex
        extends AbstractChangeable
        implements AlphabetIndex, java.io.Serializable
{
  private /*final*/ Reference alphaRef;
  private Symbol[] symbols;

  private ChangeListener indexBuilder; // do not remove these 	 
  private ChangeListener adapter;      // do not remove these
   
  // hack for bug in compaq 1.2?
  protected ChangeSupport getChangeSupport(ChangeType ct) {
    return super.getChangeSupport(ct);
  }

  public LinearAlphabetIndex(FiniteAlphabet alpha) {
    // lock the alphabet
    alpha.addChangeListener(ChangeListener.ALWAYS_VETO, Alphabet.SYMBOLS);

    this.alphaRef = new WeakReference(alpha);

    this.symbols = buildIndex(alpha);

    //these change listeners are needed to rebuild the index when the alphabet changes. 
    alpha.addChangeListener( 	 
       indexBuilder = new IndexRebuilder(), 	 
       Alphabet.SYMBOLS 	 
     ); 	 
  	 
     this.addChangeListener( 	 
       adapter = new ChangeAdapter() { 	 
         public void postChange(ChangeEvent ce) { 	 
           symbols = (Symbol[] ) ce.getChange(); 	 
         } 	 
       }, 	 
       AlphabetIndex.INDEX 	 
     ); 	 
     
    // unlock the alphabet
    alpha.removeChangeListener(ChangeListener.ALWAYS_VETO, Alphabet.SYMBOLS);
  }

  public LinearAlphabetIndex(Symbol[] syms)
  throws BioException {
    Set si = new HashSet();
    Symbol[] symbols = new Symbol[syms.length];
    for(int i = 0; i < syms.length; i++) {
      Symbol s = syms[i];
      symbols[i] = s;
      si.add(s);
    }

    FiniteAlphabet alpha = new SimpleAlphabet(si);
    this.alphaRef = new WeakReference(alpha);
    alpha.addChangeListener(ChangeListener.ALWAYS_VETO, Alphabet.SYMBOLS);
    this.symbols = symbols;
  }

  private Symbol[] buildIndex(FiniteAlphabet alpha) {
    Symbol[] symbols = new Symbol[alpha.size()];

    int i = 0;
    Iterator s = alpha.iterator();
    while(s.hasNext()) {
      symbols[i++] = (Symbol) s.next();
    }

    return symbols;
  }

  public FiniteAlphabet getAlphabet() {
    return (FiniteAlphabet) alphaRef.get();
  }

  public Symbol symbolForIndex(int i) throws IndexOutOfBoundsException {

    try {
      return symbols[i];
    } catch (IndexOutOfBoundsException e) {
      throw new IndexOutOfBoundsException("Can't find symbol for index " + i);
    }
  }

  public int indexForSymbol(Symbol s) throws IllegalSymbolException {

    for(int i = 0; i < symbols.length; i++) {
      if(s == symbols[i]) {
        return i;
      }
    }
    getAlphabet().validate(s);
    if(s instanceof AtomicSymbol) {
      throw new BioError(
        "Assertion Failure: " +
        "Symbol " + s.getName() + " was not an indexed member of the alphabet " +
        getAlphabet().getName() + " despite being in the alphabet."
      );
    } else {
      throw new IllegalSymbolException("Symbol must be atomic to be indexed.");
    }
  }

  protected class IndexRebuilder extends ChangeForwarder {
    public IndexRebuilder() {
      super(
            LinearAlphabetIndex.this,
                LinearAlphabetIndex.this.getChangeSupport(AlphabetIndex.INDEX)
      );
    }

    public ChangeEvent generateEvent(ChangeEvent ce)
    throws ChangeVetoException {
      if(ce.getType() != Alphabet.SYMBOLS) {
        return null;
      }

      // todo: can this be dropped now? MRP
      /*
      Object change = ce.getChange();
      Object previous = ce.getPrevious();

      if( (change == null) || (previous != null) ) {
        throw new ChangeVetoException(
          ce,
          "Can not update index as either a symbol is being removed, " +
          "or the alphabet has substantialy changed"
        );
      }

      if(! (change instanceof AtomicSymbol) ) {
        throw new ChangeVetoException(
          ce,
          "Can not update index as the symbol being added is not atomic"
        );
      }
      */

      return new ChangeEvent(
        getSource(), AlphabetIndex.INDEX,
        // fixme: buildIndex should be called using the proposed new alphabet
        LinearAlphabetIndex.this.buildIndex((FiniteAlphabet) ce.getSource()),
        symbols,
        ce
      );
    }
  }
}
