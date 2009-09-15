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

package org.biojavax.bio.seq;

/**
 * Holds info about base positions.
 * @author Richard Holland
 * @since 1.5
 */
public interface Position {
        
    /**
     * The empty position lies nowhere.
     */
    public static final Position EMPTY_POSITION = new SimplePosition(false,false,0);

    /**
     * A symbol representing a position that falls in between two bases,
     * eg. 2^3 falls somewhere in the gap between 2 and 3.
     */
    public static final String BETWEEN_BASES = "^";
    
    /**
     * A symbol representing a position that occupies a single base somewhere
     * in a range, eg. 5.10 falls on some base between 5 and 10.
     */
    public static final String IN_RANGE = ".";
    
    /**
     * Returns true if the position has a fuzzy start.
     * @return the fuzziness of the start.
     */
    public boolean getFuzzyStart();    
    
    /**
     * Returns true if the position has a fuzzy end.
     * @return the fuzziness of the end.
     */
    public boolean getFuzzyEnd();
    
    /**
     * Returns the beginning of the range of bases this base could lie in.
     * If this position is a single position, then start=end.
     * @return the start of this position.
     */
    public int getStart();  
    
    /**
     * Returns the end of the range of bases this base could lie in.
     * If this position is a single position, then start=end.
     * @return the end of this position.
     */
    public int getEnd();
    
    /**
     * Takes this position and returns a copy translated by 'distance' bases.
     * @param distance the distance to translate it.
     * @return the translated position.
     */
    public Position translate(int distance);
    
    /**
     * Returns the type of this position if it is not a point/single position.
     * Types are usually BETWEEN_BASES or IN_RANGE but could be any string value.
     * @return the type of this position.
     */
    public String getType();
}
