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

import java.util.Set;

import org.biojava.bio.Annotatable;
import org.biojava.utils.ChangeVetoException;


/**
 * Annotatable objects that can have rich annotations.
 * @author Richard Holland
 * @author George Waldon
 * @see RichAnnotation
 * @since 1.5
 */
public interface RichAnnotatable extends Annotatable {

    /**
      * Return the associated annotation object.
      *
      * @return a RichAnnotation object, never null
      */
    public RichAnnotation getRichAnnotation ();

    /**
     * Returns the set of notes associated with this object. Would normally
     * delegate call to internal RichAnnotation instance.
     * @return set a set of Note objects.
     * @see Note
     */
    public Set getNoteSet();
    
    /**
     * Clears the notes associated with this object and replaces them with
     * the contents of this set. Would normally delegate call to internal
     * RichAnnotation instance.
     * @param notes the set of Note objects to replace the existing ones with.
     * @throws ChangeVetoException if the set is null or contains any objects
     * that are not Note objects.
     * @see Note
     */
    public void setNoteSet(Set notes) throws ChangeVetoException;
    
}
