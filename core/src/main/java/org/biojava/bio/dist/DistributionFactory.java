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

import java.io.NotSerializableException;
import java.io.ObjectStreamException;
import java.io.Serializable;

import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.utils.StaticMemberPlaceHolder;

/**
 * <p>
 * A thing that can make Distributions.
 * </p>
 *
 * <p>
 * This decouples programs from needing to know what implementation of Distribution
 * to instantiate for a given alphabet. It also lets you parameterise model creation
 * for things like profile HMMs.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.0
 */
public interface DistributionFactory {
  /**
   * Generate a new Distribution as requested.
   *
   * @param alpha  the emission alphabet for the state
   * @return a new Distribution instance
   * @throws IllegalAlphabetException if the factory is unable to generate a
   *         distribution for the required alphabet
   */
  Distribution createDistribution(Alphabet alpha)
  throws IllegalAlphabetException;
  
  /**
   * <p>
   * The default DistributionFactory object.
   * </p>
   *
   * <p>
   * You may wish to alias this within your scripts with something like:
   * DistributionFactory dFact = DistributionFactory.DEFAULT; dFact.createDistribution(...);
   * </p>
   */
  static DistributionFactory DEFAULT = new DefaultDistributionFactory();
  
  /**
   * <p>
   * The default DistributionFactory implementation.
   * </p>
   *
   * <p>
   * It knows about hand-optimized implementations for some alphabets (like DNA)
   * without the optimized classes needing to be exposed from the DP package.
   * </p>
   *
   * @author Matthew Pocock
   */
  class DefaultDistributionFactory implements DistributionFactory, Serializable {
    public Distribution createDistribution(Alphabet alpha)
    throws IllegalAlphabetException {
      Distribution dis;
      if(! (alpha instanceof FiniteAlphabet) ) {
        throw new IllegalAlphabetException(
          "The default StateFactory implementation can only produce states over " +
          "finite alphabets, not " + alpha.getName()
        );
      }
      FiniteAlphabet fa = (FiniteAlphabet) alpha;
      
      //if(fa == DNATools.getDNA()) {
      //  dis = new DNADistribution();
      //} else {
        dis = new SimpleDistribution(fa);
      //}
    
      return dis;
    }
    
    private Object writeReplace() throws ObjectStreamException {
      try {
        return new StaticMemberPlaceHolder(DistributionFactory.class.getField("DEFAULT"));
      } catch (NoSuchFieldException nsfe) { 
        throw new NotSerializableException(nsfe.getMessage());
      }
    }
  }
}
