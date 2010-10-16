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

import java.io.NotSerializableException;
import java.io.ObjectStreamException;
import java.io.Serializable;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.seq.io.IntegerTokenization;
import org.biojava.bio.seq.io.SubIntegerTokenization;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.SingletonList;
import org.biojava.utils.StaticMemberPlaceHolder;
import org.biojava.utils.Unchangeable;
import org.biojava.utils.cache.WeakValueHashMap;

/**
 * <p>
 * An efficient implementation of an Alphabet over the infinite set of integer
 * values.
 * </p>
 *
 * <p>
 * This class can be used to represent lists of integer numbers as a
 * SymbolList with the alphabet IntegerAlphabet. These lists can then be
 * annotated with features, or fed into dynamic-programming algorithms, or
 * processed as per any other SymbolList object.
 * </p>
 *
 * <p>
 * Object identity should be used to decide if two IntegerSymbol objects are
 * the same. IntegerAlphabet ensures that all IntegerSymbol instances are
 * canonicalized.
 * </p>
 *
 * @author Matthew Pocock
 * @author Mark Schreiber
 * @author Thomas Down
 */

public final class IntegerAlphabet
  extends
    Unchangeable
  implements
    Alphabet,
    Serializable
{
  /**
   * The singleton instance of the IntegerAlphabet class.
   */
  public static IntegerAlphabet INSTANCE;

  private Object writeReplace() throws ObjectStreamException {
    try {
      return new StaticMemberPlaceHolder(IntegerAlphabet.class.getField("INSTANCE"));
    } catch (NoSuchFieldException nsfe) {
      throw new NotSerializableException(nsfe.getMessage());
    }
  }

  /**
   * Construct a finite contiguous subset of the <code>IntegerAlphabet</code>.
   * Useful for making CrossProductAlphabets with other <code>FiniteAlphabet</code>s.
   *
   * @param min the lower bound of the Alphabet
   * @param max the upper bound of the Alphabet
   * @throws IllegalArgumentException if max < min
   * @return A FiniteAlphabet from min to max <b>inclusive</b>.
   */
  public static SubIntegerAlphabet getSubAlphabet(int min, int max)
  throws IllegalArgumentException {
    String name = "SUBINTEGER["+min+".."+max+"]";
    if(AlphabetManager.registered(name)){
      return (SubIntegerAlphabet) (AlphabetManager.alphabetForName(name));
    }
          
    FiniteAlphabet a = new SubIntegerAlphabet(min, max);
    AlphabetManager.registerAlphabet(a.getName(),a);
  
    return (SubIntegerAlphabet) (AlphabetManager.alphabetForName(name));
  }

  /**
   * Retrieve a SymbolList view of an array of integers.
   * <p>
   * The returned object is a view onto the underlying array, and does not copy
   * it. Changes made to the original array will alter the symulting SymbolList.
   *
   * @param iArray  the array of integers to view
   * @return a SymbolList over the IntegerAlphabet that represent the values in
   *         iArray
   */
  public static SymbolList fromArray(int [] iArray) {
    return new IntegerArray(iArray);
  }

  /**
   * Retrieve the single IntegerAlphabet instance.
   *
   * @return the singleton IntegerAlphabet instance
   */
  public static IntegerAlphabet getInstance() {
    if(INSTANCE == null) {
      INSTANCE = new IntegerAlphabet();
      //add an alias
      AlphabetManager.registerAlphabet("Alphabet of all integers.", INSTANCE);
    }

    return INSTANCE;
  }

  /**
   * Canonicalization map for ints and references to symbols.
   */
    private WeakValueHashMap intToSym;

  private IntegerAlphabet() {
      intToSym = new WeakValueHashMap();
  }

  /**
   * Retrieve the Symbol for an int.
   *
   * @param val  the int to view
   * @return a IntegerSymbol embodying val
   */

  public synchronized IntegerSymbol getSymbol(int val) {
      Integer i = new Integer(val);
      IntegerSymbol sym = (IntegerSymbol) intToSym.get(i);
      if(sym == null) {
          sym = new IntegerSymbol(val);
          intToSym.put(i, sym);
      }
      return sym;
  }

  public Symbol getGapSymbol() {
    return AlphabetManager.getGapSymbol(getAlphabets());
  }

  public Annotation getAnnotation() {
    return Annotation.EMPTY_ANNOTATION;
  }

  public List getAlphabets() {
    return new SingletonList(this);
  }

  public Symbol getSymbol(List symList)
  throws IllegalSymbolException {
    throw new BioError("Unimplemneted method");
  }

  public Symbol getAmbiguity(Set symSet)
  throws IllegalSymbolException {
    throw new BioError("Unimplemneted method");
  }

  public boolean contains(Symbol s) {
    if(s instanceof IntegerSymbol) {
      return true;
    } else {
      return false;
    }
  }

  public void validate(Symbol s) throws IllegalSymbolException {
    if(!contains(s)) {
      throw new IllegalSymbolException(
        "Only symbols of type IntegerAlphabet.IntegerSymbol are valid for this alphabet.\n" +
        "(" + s.getClass() + ") " + s.getName()
      );
    }
  }

  public String getName() {
    return "INTEGER";
  }

  /**
   * Creates a new parser (Mark Schreiber 3 May 2001).
   *
   * @param name Currently only "token" is supported. You may also
   * use "default" as a synonym of "token"
   * @return an IntegerParser.
   */
  public SymbolTokenization getTokenization(String name) {
    if(name.equals("token") || name.equals("default")){
      return new IntegerTokenization();
    }else{
      throw new NoSuchElementException(name + " parser not supported by IntegerAlphabet yet");
    }
  }

  /**
   * A single int value.
   * <p>
   * @author Matthew Pocock
   */
  public static class IntegerSymbol
    extends
      Unchangeable
    implements
      AtomicSymbol,
      Serializable
  {
    private final int val;
    private final Alphabet matches;

    public Annotation getAnnotation() {
      return Annotation.EMPTY_ANNOTATION;
    }

    public String getName() {
      return val + "";
    }

    public int intValue() {
      return val;
    }

    public Alphabet getMatches() {
      return matches;
    }

    public List getSymbols() {
      return new SingletonList(this);
    }

    public Set getBases() {
      return Collections.singleton(this);
    }

    protected IntegerSymbol(int val) {
      this.val = val;
      this.matches = new SingletonAlphabet(this);
    }

    public int hashCode(){
      int result = 17;
      result = 37*result+intValue();
      return result;
    }

    public boolean equals(Object o){
      if(o == this) return true;
      if(o instanceof IntegerSymbol){
        IntegerSymbol i = (IntegerSymbol) o;
        if (i.intValue() == this.intValue()) {
          return true;
        }
      }
      return false;
    }
  }

  /**
   * A light-weight implementation of SymbolList that allows an array to
   * appear to be a SymbolList.
   *
   * @author Matthew Pocock
   */
  private static class IntegerArray
  extends AbstractSymbolList implements Serializable {
    private final int [] iArray;

    public Alphabet getAlphabet() {
      return INSTANCE;
    }

    public Symbol symbolAt(int i) {
      return new IntegerSymbol(iArray[i-1]);
    }

    public int length() {
      return iArray.length;
    }

    public IntegerArray(int [] iArray) {
      this.iArray = iArray;
    }
  }

  /**
   * A class to represent a finite contiguous subset of the infinite IntegerAlphabet
   *
   * @author Mark Schreiber
   * @author Matthew Pocock
   * @since 1.3
   */
  public static class SubIntegerAlphabet
  extends AbstractAlphabet {
    private int min;
    private int max;
    private String name; // cache this for performance

    /**
     * Construct a contiguous sub alphabet with the integers from min to max inclusive.
     */
    private SubIntegerAlphabet(int min, int max) throws IllegalArgumentException{
      if(max < min) {
        throw new IllegalArgumentException(
          "min must be less than max: " +
          min + " : " + max
        );
      }

      this.min = min;
      this.max = max;

      this.name = "SUBINTEGER["+min+".."+max+"]";
    }

    public String getName() {
      return name;
    }

    protected boolean containsImpl(AtomicSymbol sym) {
      if(!IntegerAlphabet.getInstance().contains(sym)) {
        return false;
      }

      IntegerSymbol is = (IntegerSymbol) sym;
      return is.intValue() >= min && is.intValue() <= max;
    }

    /**
     * @param name Currently only "token" is supported.
     * @return an IntegerParser.
     */
    public SymbolTokenization getTokenization(String name) {
      if(name.equals("token") || name.equals("default")){
        return new SubIntegerTokenization(this);
      }else{
        throw new NoSuchElementException(name + " parser not supported by IntegerAlphabet yet");
      }
    }

    public IntegerSymbol getSymbol(int val)
    throws IllegalSymbolException {
      if(val < min || val > max) {
        throw new IllegalSymbolException(
          "Could not get Symbol for value " +
          val + " as it is not in the range " +
          min + " : " + max
        );
      }

      return IntegerAlphabet.getInstance().getSymbol(val);
    }

    public int size() {
      return max - min + 1;
    }

    public List getAlphabets() {
      return new SingletonList(this);
    }


    protected AtomicSymbol getSymbolImpl(List symL) throws
        IllegalSymbolException {

      if (symL.size() != 1) {
        throw new IllegalSymbolException(
            "SubIntegerAlphabet is one-dimensional: " + this.getName() +
            " : " + symL);
      }

      AtomicSymbol s = (AtomicSymbol) symL.get(0);
      this.validate(s);
      return s;
    }

    protected void addSymbolImpl(AtomicSymbol sym)
    throws ChangeVetoException {
      throw new ChangeVetoException(
        "Can't add symbols to immutable alphabet " +
        getName()
      );
    }

    public void removeSymbol(Symbol sym)
    throws ChangeVetoException {
      throw new ChangeVetoException(
        "Can't remove symbols from immutable alphabet " +
        getName()
      );
    }

    public Iterator iterator() {
      return new Iterator() {
        int indx = min;

        public boolean hasNext() {
          return indx <= max;
        }

        public Object next() {
          try {
            Symbol sym = getSymbol(indx);
            indx++;
            return sym;
          } catch (IllegalSymbolException ise) {
            throw new BioError(
              "Assertion Failure: symbol " + indx +
              " produced by iterator but not found in " + getName()
              ,ise
            );
          }
        }

        public void remove() {
          throw new UnsupportedOperationException();
        }
      };
    }

    public Annotation getAnnotation() {
      return Annotation.EMPTY_ANNOTATION;
    }
  }
}
