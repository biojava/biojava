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

import java.util.ArrayList;
import java.util.Random;

import org.biojava.bio.BioError;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;

/**
 * <p>Simple Implementation of the <code>CrossOverFunction</code> interface</p>
 *
 * @author Mark Schreiber
 * @version 1.0
 * @since 1.5
 */

public class SimpleCrossOverFunction extends AbstractCrossOverFunction {

  public SimpleCrossOverFunction() {
  }


  // This is the one that actually does the work
  public GACrossResult performCrossOver(SymbolList chromA,
                                        SymbolList chromB)
    throws ChangeVetoException {

    ArrayList crossPoints = new ArrayList();
    Random rand = new Random();

    //do the actual crossing!
    double crossProb;
    //don't use <= crhomA.length() as there is no point crossing at the last pos
    for (int i = 1; i < chromA.length() && i < chromB.length(); i++) {

      //crossOverProbs might be shorter than i, in this case use the last prob
      if(i - 1 > getCrossOverProbs().length -1)
        crossProb = getCrossOverProbs()[getCrossOverProbs().length -1];
      else
        crossProb = getCrossOverProbs()[i-1];

      if(crossPoints.size() >= getMaxCrossOvers())
        break;

      if(rand.nextDouble() <= crossProb){


        //record the cross
        crossPoints.add(new PointLocation(i));

        //do a cross over
        SymbolList aReplace = chromB.subList(i, chromB.length());
        SymbolList bReplace = chromA.subList(i, chromA.length());

        //replace chromA from cross point down with chromB from crosspoint down
        Edit ed = new Edit(i, aReplace.length(), aReplace);
        try {
          chromA.edit(ed);
        }catch (IllegalAlphabetException ex) {
          //can't happen
          throw new BioError(ex);
        }

        //do the reciprocal
        ed = new Edit(i, bReplace.length(), bReplace);
        try {
          chromB.edit(ed);
        }
        catch (IllegalAlphabetException ex) {
          throw new BioError(ex);
        }
      }
    }


    PointLocation[] crosses = new PointLocation[crossPoints.size()];
    crosses = (PointLocation[])crossPoints.toArray(crosses);
    return new SimpleGACrossResult(crosses, new SymbolList[]{chromA, chromB});
  }
}