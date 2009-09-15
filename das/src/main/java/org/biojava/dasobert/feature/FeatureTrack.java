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
 * Created on Feb 9, 2005
 *
 */
package org.biojava.dasobert.feature;

import java.util.List;

/**
 * A feature corresponds to a track in Ensembl
 * 
 * @author Andreas Prlic
 *
 */
public interface FeatureTrack {
    
	
	public Object clone();
	
    /** returns true if the specified sequence position is within the range of this Feature
     * 
     * @param seqPosition the position to check
     * @return true if the position is within the ranges of the segments of this feature
     */
    public boolean overlaps(int seqPosition);
    
    public  String toString();

    public  void setSource(String s);

    public  String getSource();

    public  void setName(String nam);

    public  String getName();

    public  void setMethod(String methd);

    public  String getMethod();

    public  void setType(String typ);

    public  String getType();

    public  void setNote(String nte);

    public  String getNote();

    public  void setLink(String lnk);

    public  String getLink();
    
    public  void setScore(String score);
    
    public  String getScore();
    
    public void setOrientation(String orientation);
    
    public String getOrientation();
    
    /** test if two features are equivalent
     * 
     * @param feat feature to compare with 
     * @return true if equivalend
     */
    public abstract boolean equals(FeatureTrack feat);

    /** add a segment to this feature
     * 
     * @param start position
     * @param end position 
     * @param name of feature
     */
    public abstract void addSegment(int start, int end, String name);

    public abstract void addSegment(Segment s);

    public abstract List<Segment> getSegments();
    
    /** set the data from the DAS - type - id field
     * (used for Ontology support)
     * @param typeID
     */
    public void setTypeID(String typeID);
    
    /** set the data from the DAS - type - category field
     * (used for Ontology support)
     * @param typeCategory
     */
    public void setTypeCategory(String typeCategory);
    
    public String getTypeID();
    public String getTypeCategory();
   
}