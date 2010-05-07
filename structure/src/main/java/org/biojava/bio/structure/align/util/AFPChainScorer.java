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
 * Created on Apr 7, 2010
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.util;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.SVDSuperimposer;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.jama.Matrix;

public class AFPChainScorer
{

   public static double getTMScore(AFPChain align, Atom[] ca1, Atom[] ca2) throws StructureException
   {

      // Create new arrays for the subset of atoms in the alignment.
      Atom[] ca1aligned = new Atom[align.getOptLength()];
      Atom[] ca2aligned = new Atom[align.getOptLength()];
      int pos=0;
      int[] blockLens = align.getOptLen();
      int[][][] optAln = align.getOptAln();
      assert(align.getBlockNum() <= optAln.length);
      for(int block=0;block< align.getBlockNum();block++) {
         assert(blockLens[block] <= optAln[block][0].length);
         for(int i=0;i<blockLens[block];i++) {
            ca1aligned[pos] = ca1[optAln[block][0][i]];
            ca2aligned[pos] = (Atom) ca2[optAln[block][1][i]].clone();
            pos++;
         }
      }
      assert(pos == align.getOptLength());

      //Superimpose
      SVDSuperimposer svd = new SVDSuperimposer(ca1aligned, ca2aligned);
      Matrix matrix = svd.getRotation();
      Atom shift = svd.getTranslation();

      for(Atom a : ca2aligned) {
         Calc.rotate(a, matrix);
         Calc.shift(a, shift);
      }

      return SVDSuperimposer.getTMScore(ca1aligned, ca2aligned, ca1.length, ca2.length);

   }

}
