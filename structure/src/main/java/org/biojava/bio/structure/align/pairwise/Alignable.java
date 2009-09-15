/*
 *                  BioJava development code
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
 * Created on Jan 29, 2006
 *
 */
package org.biojava.bio.structure.align.pairwise;

import org.biojava.bio.structure.align.helper.AligMatEl;
import org.biojava.bio.structure.align.helper.IndexPair;

public interface Alignable {
   public int getRows();
   public int getCols();
   public AligMatEl[][] getAligMat();
   public void setAligMat(AligMatEl[][] alignmentMatrix);
   public float getGapExtCol();
   public void setGapExtCol(float penalty);
   public float getGapExtRow();
   public void setGapExtRow(float penalty);
   public float getGapOpenCol();
   public void setGapOpenCol(float penalty);
   public float getGapOpenRow();
   public void setGapOpenRow(float penalty);
   public void setScore(float score);
   public float getScore();
   public int getPathSize();
   public void setPathSize(int pathSize);
   public void setPath(IndexPair[] path);
   public IndexPair[] getPath();
}
