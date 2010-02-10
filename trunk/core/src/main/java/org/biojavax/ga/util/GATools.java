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


package org.biojavax.ga.util;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.OrderNDistribution;
import org.biojava.bio.dist.OrderNDistributionFactory;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.seq.io.CharacterTokenization;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleAlphabet;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;


/**
 * <p> Utility methods for the GA library
 * 
 * @author Mark Schreiber
 * @version 1.0
 * @since 1.5
 */

public final class GATools {

  private static SimpleAlphabet binary;
  private static AtomicSymbol zero;
  private static AtomicSymbol one;

  static{
    zero = AlphabetManager.createSymbol("zero");
    one = AlphabetManager.createSymbol("one");

    Set syms = new HashSet();
    syms.add(zero); syms.add(one);

    binary = new SimpleAlphabet(syms, "GA_Binary");
    CharacterTokenization tk = new CharacterTokenization(binary, false);
    tk.bindSymbol(zero, '0');
    tk.bindSymbol(one, '1');

    binary.putTokenization("token", tk);

    AlphabetManager.registerAlphabet(binary.getName(), binary);
  }

  /**
   * Gets a Reference to the FlyWeight GA_Binary <code>Alphabet</code>.
   * It contains the Symbols one and zero.
   * @return the finite, flyweight Binary Alphabet
   */
  public static FiniteAlphabet getBinaryAlphabet(){
    return binary;
  }

  /**
   * @return the GA_Binary symbol "one"
   */
  public static AtomicSymbol one(){
    return one;
  }

  /**
   * Creates a <code>SymbolList</code> in the GABinary <code>Alphabet</code>
   * @param binarySequence a String like "01010000101010101" with no white space
   * @return a <code>SymbolList</code> parsed from <code>binarySequence</code>
   * @throws IllegalSymbolException if a character other than 1 or 0 is found.
   */
  public static SymbolList createBinary(String binarySequence)
       throws IllegalSymbolException{

    SymbolTokenization toke = null;
    try {
      toke =
          getBinaryAlphabet().getTokenization("token");
    }
    catch (BioException ex) {
      throw new BioError("Cannot make binary tokenization", ex);
    }

    return new SimpleSymbolList(toke, binarySequence);
  }

  /**
   * @return the GA_Binary symbol "zero"
   */
  public static AtomicSymbol zero(){
    return zero;
  }

  /**
   * Makes a 1st order distribution which is infact uniform (equivalent to a
   * uniform zero order distribution).
   * @param a the zero order Alphabet which will be multiplied into the 1st order alphabet
   * @return the "1st order" distribution
   * @throws IllegalAlphabetException if the Distribution cannot be constructed from <code>a</code>.
   */
  public static OrderNDistribution uniformMutationDistribution(FiniteAlphabet a) throws IllegalAlphabetException{
    List l = Collections.nCopies(2, a);
    Alphabet alpha = AlphabetManager.getCrossProductAlphabet(l);

    OrderNDistribution d =
        (OrderNDistribution)OrderNDistributionFactory.DEFAULT.createDistribution(alpha);

    AlphabetIndex ind = AlphabetManager.getAlphabetIndex(a);
    UniformDistribution u = new UniformDistribution(a);
    for(int i = 0; i < a.size(); i++){
      try {
        d.setDistribution(ind.symbolForIndex(i), u);
      }
      catch (IllegalSymbolException ex) {
        throw new BioError(ex); //shouldn't happen
      }
    }
    return d;
  }

  /**
   * Makes a mutation <code>Distribution</code> where the probability
   * of a <code>Symbol</code> being mutated to itself is zero and the
   * probability of it being changed to any other <code>Symbol</code> in
   * the <code>Alphabet a</code> is <code>1.0 / (a.size() - 1.0)</code>
   * @param a the <code>FiniteAlphabet</code> which mutations are sampled from.
   * @return A <code>Distribution</code> suitable for use in a <code>MutationFunction</code>
   * @throws IllegalAlphabetException if the <code>Distribution</code> cannot be made
   * over the <code>FiniteAlphabet</code>
   */
  public static OrderNDistribution standardMutationDistribution(FiniteAlphabet a) throws IllegalAlphabetException{
    List l = Collections.nCopies(2, a);
    Alphabet alpha = AlphabetManager.getCrossProductAlphabet(l);

    OrderNDistribution d =
        (OrderNDistribution)OrderNDistributionFactory.DEFAULT.createDistribution(alpha);

    AlphabetIndex ind = AlphabetManager.getAlphabetIndex(a);
    for(int i = 0; i < a.size(); i++){
      try {
        Distribution sub_dist = d.getDistribution(ind.symbolForIndex(i));

        AlphabetIndex ind2 = AlphabetManager.getAlphabetIndex(a);
        for (int j = 0; j < a.size(); j++){
          if(ind.symbolForIndex(i) == ind2.symbolForIndex(j)){
            sub_dist.setWeight(ind2.symbolForIndex(j), 0.0);
          }else{
            sub_dist.setWeight(ind2.symbolForIndex(j), 1.0/ (double)(a.size() -1));
          }
        }
      }catch (IllegalSymbolException ex) {
        throw new BioError(ex); //shouldn't happen
      }catch (ChangeVetoException ex){
        throw new BioError(ex); //shouldn't happen
      }
    }
    return d;
  }
}