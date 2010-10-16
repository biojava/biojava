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

package org.biojavax.ga;

/**
 * Used by a <code>GeneticAlgorithm.run()</code> method
 * to determine when the algorithm should stop</p>
 * 
 * @author Mark Schreiber
 * @version 1.0
 * @since 1.5
 */

public interface GAStoppingCriteria {

  /**
   * Determines if an Algorithm should stop spawning new generations
   * @param ga the Algorithm to test
   * @return true if it should stop, false otherwise.
   */
  public boolean stop(GeneticAlgorithm ga);

  /**
   * Simple Implementation of GAStoppingCriteria, signals
   * a <code>GeneticAlgorithm</code> to stop after n generations</p>
   *
   * Useful for pausing and seeing what the
   * state of the algorithm is at any particular time and possibly changing
   * parameters etc before calling the run() method of the
   * <code>GeneticAlgorithm</code> again.
   * 
   * @author Mark Schreiber
   * @version 1.0
   */
  public class MaximumGeneration implements GAStoppingCriteria {
    int maxGenerations;

    /**
     * Public Constructer
     * @param maxGenerations the number of generations to stop after
     */
    public MaximumGeneration(int maxGenerations){
      this.maxGenerations = maxGenerations;
    }

    /**
     * Stops the Algorithm if the iterations are >= <code>maxGenerations</code>
     * @param ga the genetic algorithm to assess for stopping
     * @return true if the <code>ga</code> should stop
     */
    public boolean stop(GeneticAlgorithm ga){
      return(ga.getGeneration() >= maxGenerations);
    }

    public int getMaxGenerations(){
      return maxGenerations;
    }
  }

}