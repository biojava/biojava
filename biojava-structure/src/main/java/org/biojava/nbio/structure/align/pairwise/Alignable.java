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
package org.biojava.nbio.structure.align.pairwise;

import org.biojava.nbio.structure.align.helper.AligMatEl;
import org.biojava.nbio.structure.align.helper.IndexPair;

public interface Alignable {
	int getRows();
	int getCols();
	AligMatEl[][] getAligMat();
	void setAligMat(AligMatEl[][] alignmentMatrix);
	float getGapExtCol();
	void setGapExtCol(float penalty);
	float getGapExtRow();
	void setGapExtRow(float penalty);
	float getGapOpenCol();
	void setGapOpenCol(float penalty);
	float getGapOpenRow();
	void setGapOpenRow(float penalty);
	void setScore(float score);
	float getScore();
	int getPathSize();
	void setPathSize(int pathSize);
	void setPath(IndexPair[] path);
	IndexPair[] getPath();
}
