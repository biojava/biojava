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
 * Created on 22.09.2004
 * @author Andreas Prlic
 *
 */

package org.biojava.dasobert.feature ;

import java.awt.Color;


/** a class to keep track of location information for a feature */
public interface Segment  {
  
	public Object clone();
	
	public  String toString();

	public  String getNote();

	public  void setNote(String note);

	public  void setStart(int strt);

	public  int getStart();

	public  void setEnd(int ed);

	public  int getEnd();

	public  void setName(String nam);

	public  String getName();

	public  void setColor(Color col);

	public  Color getColor();

	public  void setParent(FeatureTrack f);

	public  FeatureTrack getParent();

	public  void setTxtColor(String str);

	public  String getTxtColor();

	/** returns true if the specified sequence position is within the range of this Segment
	 * 
	 * @param seqPosition the position to check
	 * @return true if seqPos >= start && seqPos <= end
	 */
	public  boolean overlaps(int seqPosition);

	/** tests if two segments are overlapping
	 * 
	 * @param segment to compare with
	 * @return true if segments overlap
	 */
	public  boolean overlaps(Segment segment);

    
}
