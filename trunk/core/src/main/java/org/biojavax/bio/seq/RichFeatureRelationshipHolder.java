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
import java.util.Set;

import org.biojava.utils.ChangeVetoException;

/**
 * Holds feature relationships.
 * @author Richard Holland
 * @since 1.5
 */
public interface RichFeatureRelationshipHolder {
    
    /**
     * Adds a relationship to this feature holder.
     * @param relationship the relationship to add.
     * @throws ChangeVetoException if the relationship is unacceptable.
     */   
    public void addFeatureRelationship(RichFeatureRelationship relationship) throws ChangeVetoException;
    
    /**
     * Removes a relationship from this feature holder.
     * @param relationship the relationship to remove.
     * @throws ChangeVetoException if it cannot be removed.
     */
    public void removeFeatureRelationship(RichFeatureRelationship relationship) throws ChangeVetoException;

    /**
     * Returns the set of relationships held in this feature holder.
     * @return a set of RichFeatureRelationship objects.
     */
    public Set getFeatureRelationshipSet();
    
    /**
     * Clears the relations from this feautre holder and replaces them with a new set.
     * @param relationships the new set of features this holder should have. The set must 
     * contain only RichFeatureRelationship objects.
     * @throws ChangeVetoException if the new set could not be installed.
     */
    public void setFeatureRelationshipSet(Set relationships) throws ChangeVetoException;
}
