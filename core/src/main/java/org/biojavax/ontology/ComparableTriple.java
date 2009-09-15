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

package org.biojavax.ontology;

import java.util.Set;

import org.biojava.ontology.AlreadyExistsException;
import org.biojava.ontology.Triple;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * Comparable triples, obviously. Allows them to have descriptors. 
 * @author Richard Holland
 * @since 1.5
 */
public interface ComparableTriple extends Triple,Comparable,Changeable {

    public static final ChangeType DESCRIPTOR = new ChangeType(
            "This triple's descriptors have changed",
            "org.biojavax.ontology.ComparableTriple",
            "DESCRIPTOR"
            );
    
    /**
     * Adds a descriptor. Must not be null.
     * @param desc the descriptor to add.
     * @throws ChangeVetoException in case of objections.
     * @throws AlreadyExistsException if the descriptor already exists.
     * @throws IllegalArgumentException if the descriptor is missing.
     */
    public void addDescriptor(ComparableTerm desc) throws AlreadyExistsException, IllegalArgumentException,ChangeVetoException;
    
    /**
     * Removes a descriptor. Must not be null.
     * @return True if it did it, false if the descriptor did not exist.
     * @param desc the descriptor to remove.
     * @throws ChangeVetoException in case of objections.
     * @throws IllegalArgumentException if the descriptor is missing.
     */
    public boolean removeDescriptor(ComparableTerm desc) throws IllegalArgumentException,ChangeVetoException;
     
    /**
     * Clears the current set of descriptors and replaces it with the content of 
     * the set passed.
     * @param descriptors the set of ComparableTerm descriptors to add.
     * @throws ChangeVetoException in case of objections.
     * @see ComparableTerm
     */
    public void setDescriptors(Set descriptors) throws ChangeVetoException;
    
    /**
     * Returns all descriptors.
     * @return a set of all ComparableTerm descriptors, possibly empty.
     * @see ComparableTerm
     */
    public Set getDescriptors();
    
}