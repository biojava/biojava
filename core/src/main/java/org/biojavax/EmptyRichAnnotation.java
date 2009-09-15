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

import java.io.NotSerializableException;
import java.io.ObjectStreamException;
import java.io.Serializable;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.StaticMemberPlaceHolder;
import org.biojava.utils.Unchangeable;

/**
 * A place holder for a RichAnnotation that prevents null having to be used
 * @author Mark Schreiber
 * @author Richard Holland
 * @since 1.5
 */
public class EmptyRichAnnotation extends Unchangeable implements RichAnnotation, Serializable {
    
    private Note[] emptyNotes = new Note[]{};
    
    /**
     * {@inheritDoc} There are no properties in the Empty RichAnnotation object.
     * Calling this will return null.
     */
    public Object getProperty(Object key) throws NoSuchElementException {
        return null;
    }
    
    /**
     * {@inheritDoc} There are no properties in the Empty RichAnnotation object.
     * Calling this will return an empty array.
     */
    public Note[] getProperties(Object key) {
        return this.emptyNotes;
    }
    
    /**
     * {@inheritDoc} There are no notes in the Empty RichAnnotation object.
     * Calling this will return null.
     */
    public Note getNote(Note note){
        return null;
    }
    
    /**
     * {@inheritDoc} You can not add properties to the Empty RichAnnotation object
     * @throws ChangeVetoException whenever you call this method.
     */
    public void setProperty(Object key, Object value)
    throws ChangeVetoException {
        throw new ChangeVetoException(
                "You can not add properties to the Empty RichAnnotation object: " +
                key + " -> " + value
                );
    }
    
    /**
     * {@inheritDoc} You can not add Notes to the Empty RichAnnotation object.
     * @throws ChangeVetoException whenever you call this method.
     */
    public void setNoteSet(Set notes) throws ChangeVetoException{
        throw new ChangeVetoException(
                "You can not add Notes to the Empty RichAnnotation object");
    }
    
    /**
     * {@inheritDoc} You can not add Notes to the Empty RichAnnotation object.
     * @throws ChangeVetoException whenever you call this method.
     */
    public void addNote(Note note) throws ChangeVetoException{
        throw new ChangeVetoException(
                "You can not add Notes to the Empty RichAnnotation object");
    }
    
    /**
     * {@inheritDoc} Does nothing as it contains nothing.
     */
    public void clear() throws ChangeVetoException{ }
    
    
    /**
     * {@inheritDoc} You cannot remove properties from the Empty RichAnnotation
     * @throws ChangeVetoException whenever you call this method.
     */
    public void removeProperty(Object key)
    throws ChangeVetoException {
        throw new ChangeVetoException(
                "You cannot remove properties from the Empty RichAnnotation (!)"
                );
    }
    
    /**
     * {@inheritDoc} You cannot remove notes from the Empty RichAnnotation
     * @throws ChangeVetoException whenever you call this method.
     */
    public void removeNote(Note note)
    throws ChangeVetoException {
        throw new ChangeVetoException(
                "You cannot remove notes from the Empty RichAnnotation (!)"
                );
    }
    
    /**
     * {@inheritDoc}
     * @return always false as there are no properties
     */
    public boolean containsProperty(Object key) {
        return false;
    }
    
    /**
     * {@inheritDoc}
     * @return always false as there are no notes
     */
    public boolean contains(Note note){
        return false;
    }
    
    /**
     * {@inheritDoc}
     * @return an empty set
     */
    public Set keys() {
        return Collections.EMPTY_SET;
    }
    
    /**
     * {@inheritDoc}
     * @return an empty set
     */
    public Set getNoteSet(){
        return Collections.EMPTY_SET;
    }
    
    /**
     * {@inheritDoc}
     * @return an new Map with no entries
     */
    public Map asMap() {
        return new HashMap();
    }
    
    // For use during serialization
    private Object writeReplace() throws ObjectStreamException {
        try {
            return new StaticMemberPlaceHolder(RichAnnotation.class.getField("EMPTY_ANNOTATION"));
        } catch (NoSuchFieldException nsfe) {
            throw new NotSerializableException(nsfe.getMessage());
        }
    }
    
    /**
     * {@inheritDoc}
     * @return the hash code of a map with no entries
     */
    public int hashCode() {
        return asMap().hashCode();
    }
    
    /**
     * {@inheritDoc}
     * @return true if and only if o is an instance of
     *  this class or a descendant.
     */
    public boolean equals(Object o) {
        return (o instanceof EmptyRichAnnotation);
    }
    
}
