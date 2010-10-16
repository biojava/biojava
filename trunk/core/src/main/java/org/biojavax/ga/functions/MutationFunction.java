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

import org.biojava.bio.dist.OrderNDistribution;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * <p>A class that mutates a <code>SymbolList</code></p>
 * @author Mark Schreiber
 * @version 1.0
 * @since 1.5
 */

public interface MutationFunction extends Changeable{
   public static final double[] DEFAULT_MUTATION_PROBS = {0.0};

   public static final ChangeType MUTATION_PROBS =
      new ChangeType("mutation probabilities",MutationFunction.class,"MUTATION_PROBS");
   public static final ChangeType MUTATION_SPECTRUM =
      new ChangeType("mutation spectrum",MutationFunction.class,"MUTATION_SPECTRUM");

   public static final MutationFunction NO_MUTATION = new NoMutation();

  /**
   * Produces a new SymbolList by mutation. Each position <em>i</em> in the SymbolList <code>seq</code>
   * is mutated with probability <code>getMutationProbs[i]</code>. The new residue is selected at random
   * from the <code>Distribution mutation</code>. The use of an array of probabilities
   * allows the modelling of mutational hotspots. Position 0 in the array corresponds to the
   * probability of the first residue of <code>seq</code> mutating.
   *
   * If the length of the array defined in <code>getMutationProbs()</code> is shorter
   * than the length of the sequence the default behaivour of implementations will
   * be to apply the last probability to each subsequence residue. A single member
   * array will mutate all bases with equal probability.
   *
   * @param seq the sequence to mutate
   * @return The mutated sequence.
   * @throws IllegalAlphabetException If the <code>mutationSpectrum Distribution</code> is not
   * emitting Symbols from the same <code>Alphabet</code> as <code>seq</code>.
   * @throws IllegalSymbolException if the <code>mutationSpectrum Distribution</code> is not
   * conditioned with the same <code>Alphabet</code> as the <code>seq Alphabet</code>.
   * @throws ChangeVetoException if <code>seq</code> is unmodifiable
   */
  public SymbolList mutate(SymbolList seq)
    throws IllegalAlphabetException, ChangeVetoException, IllegalSymbolException;

  /**
   * Set the probability of a mutation occuring at a certain position
   * Position 0 in the array corresponds to the
   * probability of the first residue of <code>seq</code> mutating.
   *
   * If the length of the array defined in <code>getMutationProbs()</code> is shorter
   * than the length of the sequence the default behaivour of implementations will
   * be to apply the last probability to each subsequence residue. A single member
   * array will mutate all bases with equal probability.
   * @param mutationProbs an array of double values representing mutation probabilities
   * @throws ChangeVetoException if a ChangeListener vetoes the change.
   */
  public void setMutationProbs(double[] mutationProbs) throws ChangeVetoException;
  public double[] getMutationProbs();

  /**
   * Sets the <code>Distribution</code> of <code>Symbols</code> that will be selected
   * from when a mutation occurs. An <code>OrderNDistribution</code> is
   * used so that you can model a situation where the
   * identity of the 'mutant' <code>Symbol</code> is dependent on the original
   * <code>Symbol</code>. The primary use is to prevent <code>Symbols</code> mutating to
   * themselves. Another use would be to model transitions and transversions.
   *
   * @param mutationSpectrum the Distribution of 'mutant' bases to choose from.
   * @throws ChangeVetoException if a ChangeListener vetoes the change.
   */
  public void setMutationSpectrum(org.biojava.bio.dist.OrderNDistribution mutationSpectrum)throws ChangeVetoException;

  /**
   *
   * @return null if the Distribution has not been set.
   */
  public org.biojava.bio.dist.OrderNDistribution getMutationSpectrum();

  /**
   *
   * <p> Place Holder class that doesn't mutate its SymbolLists</p>
   * @author Mark Schreiber
   * @version 1.0
   */
  public static final class NoMutation implements MutationFunction{
    public double[] getMutationProbs(){return DEFAULT_MUTATION_PROBS;}
    public OrderNDistribution getMutationSpectrum(){return null;}
    public SymbolList mutate(SymbolList syml){return syml;}
    public void setMutationProbs(double[] muts) throws ChangeVetoException{
      throw new ChangeVetoException("Cannot setMutationProbs for a NoMutation function");
    }
    public void setMutationSpectrum(OrderNDistribution odn) throws ChangeVetoException{
      throw new ChangeVetoException("Cannot set a Mutation spectrum for a NoMutation function");
    }

    public boolean isUnchanging(ChangeType ct){return true;}
    public void removeChangeListener(ChangeListener c){};
    public void addChangeListener(ChangeListener cl){};
    public void addChangeListener(ChangeListener cl, ChangeType ct){};
    public void removeChangeListener(ChangeListener cl, ChangeType ct){};
  }
}