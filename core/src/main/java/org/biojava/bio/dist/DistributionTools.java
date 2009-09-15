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
package org.biojava.bio.dist;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceFactory;
import org.biojava.bio.seq.impl.SimpleSequenceFactory;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.BasisSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.PackedSymbolListFactory;
import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.SimpleSymbolListFactory;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolListFactory;
import org.biojava.bio.symbol.SymbolListViews;
import org.biojava.utils.AssertionFailure;
import org.biojava.utils.ChangeVetoException;
import org.xml.sax.SAXException;


/**
 * A class to hold static methods for calculations and manipulations using
 * Distributions.
 *
 * @author Mark Schreiber
 * @author Matthew Pocock
 * @since 1.2
 */

public final class DistributionTools {

  /**
   * Overide the constructer to prevent subclassing.
   */
  private DistributionTools(){}

  /**
   * Writes a Distribution to XML that can be read with the readFromXML method.
   *
   * @param d the Distribution to write.
   * @param os where to write it to.
   * @throws IOException if writing fails
   */
  public static void writeToXML(Distribution d, OutputStream os) throws IOException{
    new XMLDistributionWriter().writeDistribution(d, os);
  }

  /**
   * Read a distribution from XML.
   *
   * @param is  an InputStream to read from
   * @return  a Distribution parameterised by the xml in is
   * @throws IOException  if is failed
   * @throws SAXException if is could not be processed as XML
   */
  public static Distribution readFromXML(InputStream is)throws IOException, SAXException{
    XMLDistributionReader writer = new XMLDistributionReader();
    return writer.parseXML(is);
  }

  /**
   * Randomizes the weights of a <code>Distribution</code>.
   *
   * @param d the <code>Distribution</code> to randomize
   * @throws ChangeVetoException if the Distribution is locked
   */
  public static void randomizeDistribution(Distribution d)
    throws ChangeVetoException{
    Random rand = new Random();
    FiniteAlphabet a = (FiniteAlphabet)d.getAlphabet();
    AlphabetIndex ind = AlphabetManager.getAlphabetIndex(a);
    DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();
    dtc.registerDistribution(d);

    for(int i = 0; i < a.size(); i++){
      try {
        dtc.addCount(d,ind.symbolForIndex(i),rand.nextDouble());
      }
      catch (IllegalSymbolException ex) {
        throw new BioError("Alphabet has Illegal Symbols!!", ex);
      }
    }

    dtc.train();
  }

  /**
   * Make a distribution from a count.
   *
   * @param c the count
   * @return a Distrubution over the same <code>FiniteAlphabet</code> as <code>c</code>
   * and trained with the counts of <code>c</code>
   */
  public static Distribution countToDistribution(Count c){
    FiniteAlphabet a  = (FiniteAlphabet)c.getAlphabet();
    Distribution d = null;
    try{
      d = DistributionFactory.DEFAULT.createDistribution(a);
      AlphabetIndex index =
          AlphabetManager.getAlphabetIndex(a);
      DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();
      dtc.registerDistribution(d);

      for(int i = 0; i < a.size(); i++){
        dtc.addCount(d, index.symbolForIndex(i),
         c.getCount((AtomicSymbol)index.symbolForIndex(i)));
      }
      dtc.train();
    } catch (IllegalAlphabetException iae) {
      throw new AssertionFailure("Assertion failure: Alphabets don't match");
    }catch(IllegalSymbolException ise){
      throw new AssertionFailure("Assertion Error: Cannot convert Count to Distribution", ise);
    } catch (ChangeVetoException cve) {
      throw new AssertionFailure("Assertion failure: distributions or counts got locked.", cve);
    }
    return d;
  }

  /**
   * Compares the emission spectra of two distributions.
   *
   * @return true if alphabets and symbol weights are equal for the two distributions.
   * @throws BioException if one or both of the Distributions are over infinite alphabets.
   * @since 1.2
   * @param a A <code>Distribution</code> with the same <code>Alphabet</code> as
   * <code>b</code>
   * @param b A <code>Distribution</code> with the same <code>Alphabet</code> as
   * <code>a</code>
   */
  public static final boolean areEmissionSpectraEqual(Distribution a, Distribution b)
    throws BioException{
      //are either of the Dists infinite
      if(a.getAlphabet() instanceof FiniteAlphabet == false
          || b.getAlphabet() instanceof FiniteAlphabet == false){
        throw new IllegalAlphabetException("Cannot compare emission spectra over infinite alphabet");
      }
      //are alphabets equal?
      if(!(a.getAlphabet().equals(b.getAlphabet()))){
        return false;
      }
      //are emissions equal?
      for(Iterator i = ((FiniteAlphabet)a.getAlphabet()).iterator();i.hasNext();){
        Symbol s = (Symbol)i.next();
        if(a.getWeight(s) != b.getWeight(s)) return false;
      }
      return true;
  }

  /**
   * Compares the emission spectra of two distribution arrays.
   *
   * @return true if alphabets and symbol weights are equal for each pair
   * of distributions. Will return false if the arrays are of unequal length.
   * @throws BioException if one of the Distributions is over an infinite
   * alphabet.
   * @since 1.3
   * @param a A <code>Distribution[]</code> consisting of <code>Distributions</code>
   * over a <code>FiniteAlphabet </code>
   * @param b A <code>Distribution[]</code> consisting of <code>Distributions</code>
   * over a <code>FiniteAlphabet </code>
   */
  public static final boolean areEmissionSpectraEqual(Distribution[] a,
                                                      Distribution[] b)
    throws BioException{
      if(a.length != b.length) return false;
      for (int i = 0; i < a.length; i++) {
        if(areEmissionSpectraEqual(a[i], b[i]) == false){
          return false;
        }
      }
      return true;
    }

  /**
   * A method to calculate the Kullback-Liebler Distance (relative entropy).
   *
   * @param logBase  - the log base for the entropy calculation. 2 is standard.
   * @param observed - the observed frequence of <code>Symbols </code>.
   * @param expected - the excpected or background frequency.
   * @return  - A HashMap mapping Symbol to <code>(Double)</code> relative entropy.
   * @since 1.2
   */
  public static final HashMap KLDistance(Distribution observed,
                                   Distribution expected,
                                   double logBase){
    Iterator alpha = ((FiniteAlphabet)observed.getAlphabet()).iterator();
    HashMap kldist = new HashMap(((FiniteAlphabet)observed.getAlphabet()).size());

    while(alpha.hasNext()){
      Symbol s = (Symbol)alpha.next();
      try{
        double obs = observed.getWeight(s);
        double exp = expected.getWeight(s);
        if(obs == 0.0){
          kldist.put(s,new Double(0.0));
        }else{
          double entropy = obs * (Math.log(obs/exp))/Math.log(logBase);
          kldist.put(s,new Double(entropy));
        }
      }catch(IllegalSymbolException ise){
        ise.printStackTrace(System.err);
      }
    }
    return kldist;
  }

  /**
   * A method to calculate the Shannon Entropy for a Distribution.
   *
   * @param logBase  - the log base for the entropy calculation. 2 is standard.
   * @param observed - the observed frequence of <code>Symbols </code>.
   * @return  - A HashMap mapping Symbol to <code>(Double)</code> entropy.
   * @since 1.2
   */
  public static final HashMap shannonEntropy(Distribution observed, double logBase){
    Iterator alpha = ((FiniteAlphabet)observed.getAlphabet()).iterator();
    HashMap entropy = new HashMap(((FiniteAlphabet)observed.getAlphabet()).size());

    while(alpha.hasNext()){
      Symbol s = (Symbol)alpha.next();
      try{
        double obs = observed.getWeight(s);
        if(obs == 0.0){
         // entropy.put(s,new Double(0.0));
        }else{
          double e = -(Math.log(obs))/Math.log(logBase);
          entropy.put(s,new Double(e));
        }
      }catch(IllegalSymbolException ise){
        ise.printStackTrace(System.err);
      }
    }
    return entropy;
  }

  /**
   * Calculates the total Entropy for a Distribution. Entropies for individual
   * <code>Symbols</code> are weighted by their probability of occurence.
   * @param observed the observed frequence of <code>Symbols </code>.
   * @return the total entropy of the <code>Distribution </code>.
   */
  public static double totalEntropy(Distribution observed){
    HashMap ent = shannonEntropy(observed, 2.0);
    double totalEntropy = 0.0;
    try{
    for(Iterator i = ent.keySet().iterator(); i.hasNext();){
      Symbol sym = (Symbol) i.next();
      totalEntropy += observed.getWeight(sym)*((Double)ent.get(sym)).doubleValue();
    }
    }
    catch(Exception e){
      e.printStackTrace(System.err);
    }

    return totalEntropy;
  }

  /**
   * Calculates the total bits of information for a distribution.
   * @param observed - the observed frequence of <code>Symbols </code>.
   * @return the total information content of the <code>Distribution </code>.
   * @since 1.2
   */
  public static final double bitsOfInformation(Distribution observed){
    double totalEntropy = totalEntropy(observed);
    int size = ((FiniteAlphabet)observed.getAlphabet()).size();

    return Math.log((double)size)/Math.log(2.0) - totalEntropy;
  }

  /**
   * Equivalent to distOverAlignment(a, false, 0.0).
   *
   * @param a  the Alignment
   * @return   an array of Distribution instances representing columns of the
   *     alignment
   * @throws IllegalAlphabetException  if the alignment alphabet is not
   *    compattible
   */
  public static Distribution[] distOverAlignment(Alignment a)
      throws IllegalAlphabetException{
    return distOverAlignment(a,false,0.0);
  }

  /**
   * Creates a joint distribution.
   *
   * @throws IllegalAlphabetException if all sequences don't use the same alphabet
   * @param a the <code>Alignment </code>to build the <code>Distribution[]</code> over.
   * @param countGaps if true gaps will be included in the distributions
   * (NOT YET IMPLEMENTED!!, CURRENTLY EITHER OPTION WILL PRODUCE THE SAME RESULT)
   * @param nullWeight the number of pseudo counts to add to each distribution
   * @param cols a list of positions in the alignment to include in the joint distribution
   * @return a <code>Distribution</code>
   * @since 1.2
   */
  public static final Distribution jointDistOverAlignment(Alignment a,
                                                 boolean countGaps,
                                                 double nullWeight,
                                                 int[] cols)
  throws IllegalAlphabetException {
        List seqs = a.getLabels();
        FiniteAlphabet alpha =
          (FiniteAlphabet)((SymbolList)a.symbolListForLabel(seqs.get(0))).getAlphabet();
        for(int i = 1; i < seqs.size();i++){
                FiniteAlphabet test = (FiniteAlphabet)((SymbolList)a.symbolListForLabel(seqs.get(i))).getAlphabet();
                if(test != alpha){
                        throw new IllegalAlphabetException("Cannot Calculate jointDistOverAlignment() for alignments with"+
                        "mixed alphabets");
                }
        }
        List a_list = new ArrayList();
        for(int i=0; i<cols.length; i++){
                a_list.add(alpha);
        }
        Distribution dist;
        DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();
        dist =
          DistributionFactory.DEFAULT.
          createDistribution(AlphabetManager.getCrossProductAlphabet(a_list));
        dtc.setNullModelWeight(nullWeight);
    try{

        dtc.registerDistribution(dist);
        Location loc= new PointLocation(cols[0]);
        for (int j = 0; j < cols.length; j++)
            {
                Location lj = new PointLocation(cols[j]);
                loc = LocationTools.union(loc, lj);
            }
            Alignment subalign = a.subAlignment(new HashSet(seqs), loc);
            Iterator s_it = subalign.symbolListIterator();
        while(s_it.hasNext()){
            SymbolList syml = (SymbolList) s_it.next();
            Symbol s= SymbolListViews.orderNSymbolList(syml,syml.length()).symbolAt(1);
            if(countGaps == false && syml.toList().contains(a.getAlphabet().getGapSymbol())){
                    //do nothing, not counting gaps
            }else{
            dtc.addCount(dist,s,1.0);// count the symbol
            }
        }
        dtc.train();
    }catch(Exception e){
      e.printStackTrace(System.err);
    }
    return dist;
}
  /**
   * Creates an array of distributions, one for each column of the alignment.
   *
   * @throws IllegalAlphabetException if all sequences don't use the same alphabet
   * @param a the <code>Alignment </code>to build the <code>Distribution[]</code> over.
   * @param countGaps if true gaps will be included in the distributions
   * @param nullWeight the number of pseudo counts to add to each distribution,
   * pseudo counts will not affect gaps, no gaps, no gap counts.
   * @return a <code>Distribution[]</code> where each member of the array is a
   * <code>Distribution </code>of the <code>Symbols </code>found at that position
   * of the <code>Alignment </code>.
   * @since 1.2
   */
  public static final Distribution[] distOverAlignment(Alignment a,
                                                 boolean countGaps,
                                                 double nullWeight)
  throws IllegalAlphabetException {

    List seqs = a.getLabels();

    FiniteAlphabet alpha = (FiniteAlphabet)((SymbolList)a.symbolListForLabel(seqs.get(0))).getAlphabet();
    for(int i = 1; i < seqs.size();i++){
        FiniteAlphabet test = (FiniteAlphabet)((SymbolList)a.symbolListForLabel(seqs.get(i))).getAlphabet();
        if(test != alpha){
          throw new IllegalAlphabetException("Cannot Calculate distOverAlignment() for alignments with"+
          "mixed alphabets");
        }
    }

    Distribution[] pos = new Distribution[a.length()];
    DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();
    dtc.setNullModelWeight(nullWeight);

    double[] adjRatios = null;
    if(countGaps){
      adjRatios = new double[a.length()];
    }

    try{
      for(int i = 0; i < a.length(); i++){// For each position
        double gapCount = 0.0;
        double totalCount = 0.0;

        pos[i] = DistributionFactory.DEFAULT.createDistribution(alpha);
        dtc.registerDistribution(pos[i]);

        for(Iterator j = seqs.iterator(); j.hasNext();){// of each sequence
          Object seqLabel = j.next();
          Symbol s = a.symbolAt(seqLabel,i + 1);

          /*If this is working over a flexible alignment there is a possibility
          that s could be null if this Sequence is not really preset in this
          region of the Alignment. In this case it will be skipped*/
          if(s == null)
            continue;

          Symbol gap = alpha.getGapSymbol();
          if(countGaps &&
             s.equals(gap)){
             gapCount++; totalCount++;
          }else{
            dtc.addCount(pos[i],s,1.0);// count the symbol
            totalCount++;
          }
        }

        if(countGaps){
          adjRatios[i] = 1.0 - (gapCount / totalCount);
        }
      }

      dtc.train();

      if(countGaps){//need to adjust counts for gaps
        for (int i = 0; i < adjRatios.length; i++) {
          Distribution d = pos[i];
          for (Iterator iter = ((FiniteAlphabet)d.getAlphabet()).iterator();
                            iter.hasNext(); ) {
            Symbol sym = (Symbol)iter.next();
            d.setWeight(sym, (d.getWeight(sym) * adjRatios[i]));
          }
        }
      }

    }catch(Exception e){
      e.printStackTrace(System.err);
    }
    return pos;
  }


  /**
   * Creates an array of distributions, one for each column of the alignment.
   * No pseudo counts are used.
   * @param countGaps if true gaps will be included in the distributions
   * @param a the <code>Alignment </code>to build the <code>Distribution[]</code> over.
   * @throws IllegalAlphabetException if the alignment is not composed from sequences all
   *         with the same alphabet
   * @return a <code>Distribution[]</code> where each member of the array is a
   * <code>Distribution </code>of the <code>Symbols </code>found at that position
   * of the <code>Alignment </code>.
   * @since 1.2
   */
  public static final Distribution[] distOverAlignment(Alignment a,
                                                 boolean countGaps)
  throws IllegalAlphabetException {
    return distOverAlignment(a,countGaps,0.0);
  }

  /**
   * Averages two or more distributions. NOTE the current implementation ignore the null model.
   * @since 1.2
   * @param dists the <code>Distributions </code>to average
   * @return a <code>Distribution </code>were the weight of each <code>Symbol </code>
   * is the average of the weights of that <code>Symbol </code>in each <code>Distribution </code>.
   */
  public static final Distribution average (Distribution [] dists){

    Alphabet alpha = dists[0].getAlphabet();
    //check if all alphabets are the same
    for (int i = 1; i < dists.length; i++) {
      if(!(dists[i].getAlphabet().equals(alpha))){
        throw new IllegalArgumentException("All alphabets must be the same");
      }
    }

    try{
      Distribution average = DistributionFactory.DEFAULT.createDistribution(alpha);
      DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();
      dtc.registerDistribution(average);

      for (int i = 0; i < dists.length; i++) {// for each distribution
        for(Iterator iter = ((FiniteAlphabet)dists[i].getAlphabet()).iterator(); iter.hasNext(); ){//for each symbol
          Symbol sym = (Symbol)iter.next();
          dtc.addCount(average,sym,dists[i].getWeight(sym));
        }
      }


      dtc.train();
      return average;
    } catch(IllegalAlphabetException iae){//The following throw unchecked exceptions as they shouldn't happen
       throw new AssertionFailure("Distribution contains an illegal alphabet", iae);
    } catch(IllegalSymbolException ise){
       throw new AssertionFailure("Distribution contains an illegal symbol", ise);
    } catch(ChangeVetoException cve){
       throw new AssertionFailure("The Distribution has become locked", cve);
    }
  }

  /**
   * Produces a sequence by randomly sampling the Distribution.
   *
   * @param name the name for the sequence
   * @param d the distribution to sample. If this distribution is of order N a
   * seed sequence is generated allowed to 'burn in' for 1000 iterations and used
   * to produce a sequence over the conditioned alphabet.
   * @param length the number of symbols in the sequence.
   * @return a Sequence with name and urn = to name and an Empty Annotation.
   */
  public static final Sequence generateSequence(String name, Distribution d, int length){
    SymbolList sl = generateSymbolList(d, length);
    SequenceFactory fact = new SimpleSequenceFactory();
    return fact.createSequence(sl, name, name, Annotation.EMPTY_ANNOTATION);
    //return new SimpleSequence(sl,name,name,Annotation.EMPTY_ANNOTATION);
  }

/**
 * Produces a <code>SymbolList</code> by randomly sampling a Distribution.
 *
 * @param d the distribution to sample. If this distribution is of order N a
 * seed sequence is generated allowed to 'burn in' for 1000 iterations and used
 * to produce a sequence over the conditioned alphabet.
 * @param length the number of symbols in the sequence.
 * @return a SymbolList or length <code>length</code>
 */
  public static final SymbolList generateSymbolList(Distribution d, int length){
    if(d instanceof OrderNDistribution)
      return generateOrderNSymbolList((OrderNDistribution)d, length);

    SymbolList sl = null;

    List l = new ArrayList(length);
    for (int i = 0; i < length; i++) {
      l.add(d.sampleSymbol());
    }

    try {
      SymbolListFactory fact;
      if(length < 10000){
        fact = new SimpleSymbolListFactory();
      }else{
        fact = new PackedSymbolListFactory();
      }

      Symbol[] syms = new Symbol[length];
      l.toArray(syms);

      sl = fact.makeSymbolList(syms, length, d.getAlphabet());
      //sl = new SimpleSymbolList(d.getAlphabet(),l);
    }
    catch (IllegalAlphabetException ex) {
      //shouldn't happen but...
      throw new BioError("Distribution emitting Symbols not from its Alphabet?");
    }

    return sl;
  }

  private static final SymbolList generateOrderNSymbolList(OrderNDistribution d, int length){
    SymbolList sl = null;
    List l = new ArrayList(length);

    /*
     * When emitting an orderN sequence a seed sequence is required that is of the
     * length of the conditioning alphabet. The emissions will also be allowed
     * to 'burn in' for 1000 emissions so that the 'end effect' of the seed
     * is negated.
     */
     FiniteAlphabet cond = (FiniteAlphabet)d.getConditioningAlphabet();
     UniformDistribution uni = new UniformDistribution(cond);
     BasisSymbol seed = (BasisSymbol)uni.sampleSymbol();
     //using the linked list the seed becomes like a history buffer.
     LinkedList ll = new LinkedList(seed.getSymbols());

    try {

      for(int i = 0; i < 1000+ length; i++){
         //get a symbol using the seed
         Symbol sym = d.getDistribution(seed).sampleSymbol();
         if(i >= 1000){
           l.add(sym);
         }
         //add the symbol to the end of the seed
         ll.addLast(sym);
         //remove the first basis symbol of the seed
         ll.removeFirst();
         //regenerate the seed
         seed = (BasisSymbol)cond.getSymbol(ll);
       }

       SymbolListFactory fact;
       if(length < 10000){
         fact = new SimpleSymbolListFactory();
       }else{
         fact = new PackedSymbolListFactory();
       }

       Symbol[] syms = new Symbol[l.size()];
       l.toArray(syms);
       sl = fact.makeSymbolList(syms, length, d.getConditionedAlphabet());
       //sl = new SimpleSymbolList(d.getConditionedAlphabet(),l);
    }
    catch (IllegalSymbolException ex) {
      //shouldn't happen but...
      throw new BioError("Distribution emitting Symbols not from its Alphabet?",ex);
    }catch(IllegalAlphabetException ex){
      //shouldn't happen but...
      throw new BioError("Distribution emitting Symbols not from its Alphabet?",ex);
    }

    return sl;
  }

  /**
   * Generate a sequence by sampling a distribution.
   *
   * @deprecated use generateSequence() or generateSymbolList() instead.
   * @param name    the name of the sequence
   * @param d       the distribution to sample
   * @param length  the length of the sequence
   * @return        a new sequence with the required composition
   */
  protected static final Sequence generateOrderNSequence(String name, OrderNDistribution d, int length){

    SymbolList sl = generateOrderNSymbolList(d, length);
    SequenceFactory fact = new SimpleSequenceFactory();

    return fact.createSequence(sl, name, name, Annotation.EMPTY_ANNOTATION);
    //return new SimpleSequence(sl, name, name, Annotation.EMPTY_ANNOTATION);
  }

}//End of class
