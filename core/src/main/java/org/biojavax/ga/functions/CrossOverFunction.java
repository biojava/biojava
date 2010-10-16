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


package org.biojavax.ga.functions;

import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;


/**
 * Crosses two chromosomes. The basic usage of the class would be
 * something like choosing two chromosomes that you want to cross over and setting
 * these with the <code>setChromosomePair</code> method. Next you would call one of
 * the <code>performCrossOver</code> methods to do the crossing and finally you
 * would retreive the chromsome pair with the <code>getChromosomes</code> method.
 *
 * @author Mark Schreiber
 * @version 1.0
 * @since 1.5
 */

public interface CrossOverFunction extends Changeable{

  /**
   * Performs a cross between the pair of chromosomes
   * @param chromA The first chromosome in the cross
   * @param chromB The second chromosome in the cross
   * @return A <code>GACross</code> that holds the results of the cross
   * @throws ChangeVetoException if the chromosomes are unmodifiable
   */
  public GACrossResult performCrossOver(SymbolList chromA, SymbolList chromB)
    throws ChangeVetoException;

  /**
   * Sets an upper limit on the number of crosses. Its up to
   * the implementation to decide what to do when the limit is reached although
   * a good convention would be to keep only the first N crosses from the left
   * end (5' end) of the sequence.
   *
   * By convention the default upper limit is DEFAULT_MAX_CROSS (eg infinite).
   *  This value should be used as the default by all implementations.
   *
   * @param maxCrossOvers the limit on crosses
   * @throws ChangeVetoException if a ChangeListener vetoes this change
   */
  public void setMaxCrossOvers(int maxCrossOvers) throws ChangeVetoException;

  /**
   * @return the limit on crosses.
   */
  public int getMaxCrossOvers();

  /**
   * Sets the probability of crossing at each base. Each position
   * in the array corresponds to a position in the sequences to be crossed.
   *
   * The probability of a cross occuring <em>after</em> position 1 in the <code>SymbolList</code>
   * is given by <code>crossOverProbs[1]</code>. <code>CrossOverProbs[0]</code> is effectively
   * redundant as the cross would occur before the 1st position in the <code>SymbolList</code>.
   *
   * By convention if the array is shorter than the SymbolList it is being applied
   * to then the last value in the array will be applied to every subsequent residue.
   *
   * The default value in all implementations should be <code>DEFAULT_CROSS_PROB</code>
   *
   * @param crossOverProbs an array of doubles giving the probability of a
   * cross occuring at any place.
   *
   * @exception if a ChangeListener vetoes the change
   */
  public void setCrossOverProbs(double[] crossOverProbs) throws ChangeVetoException;
  public double[] getCrossOverProbs();

  public static final int DEFAULT_MAX_CROSS = Integer.MAX_VALUE;
  public static final double[] DEFAULT_CROSS_PROB = {0.0};

  public static final ChangeType MAX_CROSSES =
      new ChangeType("maximum number of crosses",CrossOverFunction.class,"MAX_CROSSES");

  public static final ChangeType CROSS_PROB =
      new ChangeType("cross over probabilities",CrossOverFunction.class,"CROSS_PROB");

  public static final CrossOverFunction NO_CROSS = new NoCross();



  /**
   * <p>A place holder CrossOverFunction that doesn't perform cross overs </p>
   * @author Mark Schreiber
   * @version 1.0
   */
  public final class NoCross implements CrossOverFunction {


    /**
     * @return a single member array equal to {0.0}
     */
    public double[] getCrossOverProbs(){
      return new double[]{0.0};
    }

    /**
     * @return 0 (after all you can't cross over with this function)
     */
    public int getMaxCrossOvers(){
        return 0;
    }

    public GACrossResult performCrossOver(SymbolList chromA, SymbolList chromB){
      return new SimpleGACrossResult(
          new PointLocation[]{},
          new SymbolList[]{chromA, chromB}
          );
    }

    public void setCrossOverProbs(double[] crossOverProb) throws ChangeVetoException{
      throw new ChangeVetoException("Cannot set the crossOverProb for a NO_CROSS function");
    }

    public void setMaxCrossOvers(int max) throws ChangeVetoException{
      throw new ChangeVetoException("Cannot change the maximum crossovers in a NO_CROSS function");
    }

    public boolean isUnchanging(ChangeType t){
      return true;
    }

    public void removeChangeListener(ChangeListener c){};
    public void addChangeListener(ChangeListener cl){};
    public void addChangeListener(ChangeListener cl, ChangeType ct){};
    public void removeChangeListener(ChangeListener cl, ChangeType ct){};
  }
}