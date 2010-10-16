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

package org.biojavax;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;
import org.biojavax.ontology.ComparableTerm;

/**
 * Note is a generic class intended to hold a term describing the note, 
 * a value to associate with that term, and a rank. It is a generic representation
 * of the various qualifier_value tables in BioSQL. It is used inside RichAnnotation
 * objects to represent annotations.
 * @author Richard Holland
 * @see RichAnnotation
 * @since 1.5
 */
public interface Note extends Comparable,Changeable {
    
    public static final ChangeType TERM = new ChangeType(
            "This note's term has changed",
            "org.biojavax.Note",
            "TERM"
            );
    public static final ChangeType RANK = new ChangeType(
            "This note's rank has changed",
            "org.biojavax.Note",
            "RANK"
            );
    public static final ChangeType VALUE = new ChangeType(
            "This note's value has changed",
            "org.biojavax.Note",
            "VALUE"
            );
    
    /**
     * Gets the term that defines this note.
     * @return a ComparableTerm object that is the key to this note.
     */
    public ComparableTerm getTerm();
    
    /**
     * Sets the term for this note. It cannot be null.
     * @param term the term to use.
     * @throws ChangeVetoException if it doesn't like the term.
     */
    public void setTerm(ComparableTerm term) throws ChangeVetoException;
    
    /**
     * Gets the value that defines this note.
     * @return a String object that is the value to this note.
     */
    public String getValue();
    
    /**
     * Sets the value for this note, or null for no value.
     * @param value the value to use.
     * @throws ChangeVetoException if it doesn't like the value.
     */
    public void setValue(String value) throws ChangeVetoException;
    
    /**
     * Gets the rank that defines this note.
     * @return an int that is the rank to this note.
     */
    public int getRank();
    
    /**
     * Sets the rank for this note.
     * @param value the rank to use.
     * @throws ChangeVetoException if it doesn't like the rank.
     */
    public void setRank(int value) throws ChangeVetoException;
    
}
