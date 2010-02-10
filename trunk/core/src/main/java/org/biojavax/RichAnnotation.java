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
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.utils.ChangeVetoException;

/**
 * An annotation collection which stores annotations as Note objects.
 * @author Richard Holland
 * @author Mark Schreiber
 * @see Note
 * @see RichAnnotatable
 * @since 1.5
 */
public interface RichAnnotation extends Annotation {
        
    public static final RichAnnotation EMPTY_ANNOTATION = new EmptyRichAnnotation();
    
    /**
     * Removes all notes from this annotation object.
     * @throws ChangeVetoException if it couldn't do it.
     */
    public void clear() throws ChangeVetoException;
    
    /**
     * Adds a note to this annotation. Must not be null.
     * If the note is already in the annotation, nothing happens.
     * @param note note to add
     * @throws ChangeVetoException if it doesn't like this.
     */
    public void addNote(Note note) throws ChangeVetoException;
    
    /**
     * Removes a note from this annotation. Must not be null.
     * If the note wasn't in the annotation, nothing happens.
     * @param note note to remove
     * @throws ChangeVetoException if it doesn't like this.
     */
    public void removeNote(Note note) throws ChangeVetoException;
    
    /**
     * Uses the term and rank to lookup a note in this annotation.
     * @param note note to lookup, using term and rank.
     * @return the matching note.
     * @throws ChangeVetoException if it doesn't like this.
     * @throws NoSuchElementException if it couldn't be found.
     */
    public Note getNote(Note note) throws NoSuchElementException;
    
    /**
     * Returns true if the given note exists in this annotation.
     * The lookup is done using the term and rank of the note.
     * @param note note to lookup
     * @return true if it is in this annotation, false if not.
     */
    public boolean contains(Note note);
    
    /**
     * Returns an immutable set of all notes in this annotation.
     * @return a set of notes.
     * @see Note
     */
    public Set getNoteSet();
    
    /**
     * Clears the notes from this annotation and replaces them with
     * those from the given set. The set cannot be null.
     * @param notes a set of Note objects to use from now on.
     * @throws ChangeVetoException if it doesn't like any of them.
     * @see Note
     */
    public void setNoteSet(Set notes) throws ChangeVetoException;
    
    /**
     * Find all the <code>Note</code>s with any rank that match the key.
     * @param key either a <code>String</code> identifier of a term from the
     * default onltology or a <code>ComparableTerm</code>
     * @return an array of matching <code>Notes</code> in order of rank or an
     * empty array if there are no matches. No implementation should ever 
     * return null!
     * @see Note
     * @see org.biojavax.ontology.ComparableTerm
     */
    public Note[] getProperties(Object key);
    
}
