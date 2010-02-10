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
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.DoubleTokenization;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.SingletonList;
import org.biojava.utils.StaticMemberPlaceHolder;
import org.biojava.utils.Unchangeable;
import org.biojava.utils.cache.WeakValueHashMap;

/**
 * <p>
 * An efficient implementation of an Alphabet over the infinite set of double
 * values.
 * </p>
 *
 * <p>
 * This class can be used to represent lists of floating-point numbers as a
 * SymbolList with the alphabet DoubleAlphabet. These lists can then be
 * annotated with features, or fed into dynamic-programming algorithms, or
 * processed as per any other SymbolList object.
 * </p>
 *
 * <p>
 * Object identity should be used to decide if two DoubleResidue objects are
 * the same. DoubleAlphabet ensures that all DoubleAlphabet instances are
 * canonicalized.
 * </p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */

public final class DoubleAlphabet
  extends
    Unchangeable
  implements
    Alphabet,
    Serializable
{

  public static DoubleAlphabet INSTANCE;

  private Object writeReplace() throws ObjectStreamException {
    try {
      return new StaticMemberPlaceHolder(DoubleAlphabet.class.getField("INSTANCE"));
    } catch (NoSuchFieldException nsfe) {
      throw new NotSerializableException(nsfe.getMessage());
    }
  }

  /**
   * <p>
   * Retrieve a SymbolList view of an array of doubles.
   * </p>
   *
   * <p>
   * The returned object is a view onto the underlying array, and does not copy
   * it. Changes made to the original array will alter the resulting SymbolList.
   * </p>
   *
   * @param dArray  the array of doubles to view
   * @return a SymbolList over the DoubleAlphabet that represent the values in
   *         dArray
   */
  public static SymbolList fromArray(double [] dArray) {
    return new DoubleArray(dArray);
  }

  /**
   * Retrieve the single DoubleAlphabet instance.
   *
   * @return the singleton DoubleAlphabet instance
   */
  public static DoubleAlphabet getInstance() {
    if(INSTANCE == null) {
      INSTANCE = new DoubleAlphabet();
      //add an alias
      AlphabetManager.registerAlphabet("Alphabet of all doubles.", INSTANCE);
    }

    return INSTANCE;
  }

    private List alphabets = null;
    private WeakValueHashMap doubleToSym;

    private DoubleAlphabet() {
        doubleToSym = new WeakValueHashMap();
    }

  /**
   * Retrieve the Symbol for a double.
   *
   * @param val  the double to view
   * @return a DoubleSymbol embodying val
   */
  public DoubleSymbol getSymbol(double val) {
      Double d = new Double(val);
      DoubleSymbol sym = (DoubleSymbol) doubleToSym.get(d);
      if (sym== null) {
          sym = new DoubleSymbol(val);
          doubleToSym.put(d, sym);
      }
      return sym;
  }

  /**
   * Retrieve the symbol for a range of doubles.
   *
   * @param minVal  the minimum value
   * @param maxVal  that maximum value
   * @return a DoubleRange containing all doubles between min and max value.
   */
  public DoubleRange getSymbol(double minVal, double maxVal) {
    // fixme: we should probably fly-weight these
    return new DoubleRange(minVal, maxVal);
  }

  public static SubDoubleAlphabet getSubAlphabet(double min, double max) {
    String name = "SUBDOUBLE["+ min +".."+ max +"]";
    if(! AlphabetManager.registered(name)){
      AlphabetManager.registerAlphabet(name, new SubDoubleAlphabet(min, max));
    }
    return (SubDoubleAlphabet)AlphabetManager.alphabetForName(name);
  }

  public Annotation getAnnotation() {
    return Annotation.EMPTY_ANNOTATION;
  }

  public boolean contains(Symbol s) {
    if(s instanceof DoubleSymbol) {
      return true;
    } else {
      return false;
    }
  }

  public void validate(Symbol s) throws IllegalSymbolException {
    if(!contains(s)) {
      throw new IllegalSymbolException(
        "Only symbols of type DoubleAlphabet.DoubleSymbol are valid for this alphabet.\n" +
        "(" + s.getClass() + ") " + s.getName()
      );
    }
  }

  public List getAlphabets() {
    if(alphabets == null) {
      alphabets = new SingletonList(this);
    }
    return alphabets;
  }

  public Symbol getGapSymbol() {
    return AlphabetManager.getGapSymbol(getAlphabets());
  }

  public Symbol getAmbiguity(Set syms) throws IllegalSymbolException {
    for(Iterator i = syms.iterator(); i.hasNext(); ) {
      Symbol sym = (Symbol) i.next();
      validate(sym);
    }
    throw new BioError("Operation not implemented");
  }

  public Symbol getSymbol(List symList) throws IllegalSymbolException {
    if(symList.size() != 1) {
      throw new IllegalSymbolException(
        "Can't build symbol from list " + symList.size() + " long"
      );
    }

    Symbol s = (Symbol) symList.get(0);
    validate(s);
    return s;
  }

  public String getName() {
    return "DOUBLE";
  }

  public SymbolTokenization getTokenization(String name) {
    if(!name.equals("name")) {
        throw new NoSuchElementException(
          "No parsers supported by DoubleAlphabet called " + name
        );
    }
    return new DoubleTokenization();
  }

  /**
   * A single double value.
   *
   *  Get these via <code>DoubleAlphabet.getSymbol(double)</code>.
   *
   * @author Matthew Pocock
   */
  public static class DoubleSymbol
    extends
      Unchangeable
    implements
      AtomicSymbol,
      Serializable
  {
    private final double val;
    private final Alphabet matches;

    public Annotation getAnnotation() {
      return Annotation.EMPTY_ANNOTATION;
    }

    public String getName() {
      return val + "";
    }

    /**
     * @return the double value associated with this double symbol
     */
    public double doubleValue() {
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

    private DoubleSymbol(double val) {
      this.val = val;
      this.matches = new SingletonAlphabet(this);
    }
  }

  /**
   * A range of double values.
   *
   *  Get these via <code>DoubleAlphabet.getSymbol(double, double)</code>.
   *
   * @author Matthew Pocock
   */
  public static class DoubleRange
    extends
      Unchangeable
    implements
      BasisSymbol,
      Serializable
  {
    private final double minVal;
    private final double maxVal;
    private final Alphabet matches;

    public Annotation getAnnotation() {
      return Annotation.EMPTY_ANNOTATION;
    }

    public String getName() {
      return "DoubleRange["+ minVal +".."+ maxVal +"]";
    }

    public Alphabet getMatches() {
      return matches;
    }

    public List getSymbols() {
      return Arrays.asList(new Symbol[] { this });
    }

    public double getMinValue() {
      return minVal;
    }

    public double getMaxValue() {
      return maxVal;
    }

    protected DoubleRange(double minVal, double maxVal) {
      this.minVal = minVal;
      this.maxVal = maxVal;
      this.matches = DoubleAlphabet.getSubAlphabet(minVal, maxVal);
    }
  }

  /**
   * A light-weight implementation of SymbolList that allows an array to
   * appear to be a SymbolList.
   *
   * @author Matthew Pocock
   */
  private static class DoubleArray
  extends
    AbstractSymbolList
  implements
    Serializable
  {
    private final double [] dArray;

    public Alphabet getAlphabet() {
      return INSTANCE;
    }

    public Symbol symbolAt(int i) {
      return new DoubleSymbol(dArray[i-1]);
    }

    public int length() {
      return dArray.length;
    }

    public DoubleArray(double [] dArray) {
      this.dArray = dArray;
    }
  }

  /**
   * A class to represent a contiguous range of double symbols.
   *
   * @author Matthew Pocock
   */
  public static class SubDoubleAlphabet
  extends
    Unchangeable
  implements
    Alphabet, Serializable
  {
    private final double min;
    private final double max;
    private final String name;

    private SubDoubleAlphabet(double min, double max) {
      this.min = min;
      this.max = max;
      this.name = "SUBDOUBLE["+ min +".."+ max +"]";
    }

    /**
     * To prevent duplication of a what should be a
     * single instance of an existing alphabet. This method
     * was written as protected so that subclasses even from
     * other packages will inherit it. It should only be overridden
     * with care.
     */
    protected Object readResolve() throws ObjectStreamException {
      try {
        return AlphabetManager.alphabetForName(this.getName());
      }
      catch (NoSuchElementException nse) {
        //a custom alphabet has been sent to your VM, register it.
        AlphabetManager.registerAlphabet(this.getName(), this);
        return this;
      }
    }


    public String getName() {
      return name;
    }

    public Annotation getAnnotation() {
      return Annotation.EMPTY_ANNOTATION;
    }

    public List getAlphabets() {
      return Arrays.asList(new Alphabet[] { this });
    }

    public Symbol getSymbol(List rl)
    throws IllegalSymbolException {
      if(rl.size() != 1) {
        throw new IllegalSymbolException(
          "SubDoubleAlphabet is one-dimensional: " + this.getName() +
          " : " + rl );
      }

      Symbol s = (Symbol) rl.get(0);

      validate(s);

      return s;
    }

    public DoubleSymbol getSymbol(double val) throws IllegalSymbolException {
      if (val < min || val > max) {
        throw new IllegalSymbolException(
            "Could not get Symbol for value " +
            val + " as it is not in the range " +
            min + " : " + max
            );
      }

      return DoubleAlphabet.getInstance().getSymbol(val);
}


    public Symbol getAmbiguity(Set syms) {
      throw new BioError("Operation not implemented");
    }

    public Symbol getGapSymbol() {
      return getInstance().getGapSymbol();
    }

    public boolean contains(Symbol s) {
      if(s instanceof DoubleSymbol) {
        double val = ((DoubleSymbol) s).doubleValue();
        if(val >= min && val <= max) {
          return true;
        }
      }

      if(s instanceof DoubleRange) {
        DoubleRange dr = (DoubleRange) s;
        if(dr.getMinValue() >= min || dr.getMaxValue() <= max) {
          return true;
        }
      }

      return false;
    }

    public void validate(Symbol sym)
    throws IllegalSymbolException {
      if(!contains(sym)) {
        throw new IllegalSymbolException(
          "This alphabet " + this.getName() +
          " does not contain the symbol " + sym );
      }
    }

    public SymbolTokenization getTokenization(String name)
    throws BioException {
      return getInstance().getTokenization(name);
    }
  }
}
