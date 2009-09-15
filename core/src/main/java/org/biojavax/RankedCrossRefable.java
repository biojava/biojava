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

import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * Defines an object as being able to have ranked cross references associated
 * with it. 
 * @author Richard Holland
 * @see RankedCrossRef
 * @since 1.5
 */
public interface RankedCrossRefable extends Changeable {
    
    /**
     * Returns the set of all ranked cross references associated with an object.
     * @return a set of RankedCrossRef objects.
     */
    public Set getRankedCrossRefs();
    
    /** 
     * Sets the ranked cross references associated with an object. Null will 
     * throw an exception but the empty set is fine.
     * @param crossrefs a set of RankedCrossRef objects.
     * @throws ChangeVetoException if the set was null or any of its contents
     * were not RankedCrossRef objects.
     */
    public void setRankedCrossRefs(Set crossrefs) throws ChangeVetoException;
    
    /**
     * Adds a ranked cross reference to the existing set. If already present, this
     * call is ignored. Null values are not acceptable.
     * @param crossref the ranked cross reference to add.
     * @throws ChangeVetoException if the parameter is null.
     */
    public void addRankedCrossRef(RankedCrossRef crossref) throws ChangeVetoException;
    
    /**
     * Removes a ranked cross reference from the existing set. If not present, this
     * call is ignored. Null values are not acceptable.
     * @param crossref the ranked cross reference to remove.
     * @throws ChangeVetoException if the parameter is null.
     */
    public void removeRankedCrossRef(RankedCrossRef crossref) throws ChangeVetoException;
}
