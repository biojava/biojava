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

import java.util.Random;

import org.biojava.bio.dist.OrderNDistribution;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;

/**
 * <p> Simple no frills Implementation of the MutationFunction interface</p>
 * <p> This class is final, custom implementations should extend <code>
 * AbstractMutationFunction</code></p>
 *
 * @author Mark Schreiber
 * @version 1.0
 * @since 1.5
 */

public final class SimpleMutationFunction extends AbstractMutationFunction {

  public SimpleMutationFunction() {
  }

  public SymbolList mutate(SymbolList seq)
      throws ChangeVetoException, IllegalAlphabetException, IllegalSymbolException {

    int maxIndex = getMutationProbs().length -1;
    OrderNDistribution d = getMutationSpectrum();
    Random r = new Random();

    for (int i = 1; i < seq.length(); i++) {
      int index = Math.min(i-1, maxIndex);
      double mutProb = getMutationProbs()[index];

      if(r.nextDouble() < mutProb){

        Edit e = new Edit(i, seq.getAlphabet(),
                          d.getDistribution(seq.symbolAt(i)).sampleSymbol());
        seq.edit(e);

      }
    }

    return seq;
  }

}