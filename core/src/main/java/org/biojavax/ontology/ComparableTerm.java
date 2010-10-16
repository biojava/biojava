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

import org.biojava.ontology.Term;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;
import org.biojavax.RankedCrossRefable;

/**
 * Makes Term objects comparable properly and adds some extra features to them.
 * @author Richard Holland
 * @since 1.5
 */
public interface ComparableTerm extends Term,RankedCrossRefable,Comparable,Changeable {
    
    public static final ChangeType IDENTIFIER = new ChangeType(
            "This term's identifier has changed",
            "org.biojavax.ontology.ComparableTerm",
            "IDENTIFIER"
            );
    public static final ChangeType OBSOLETE = new ChangeType(
            "This term's obsolescence has changed",
            "org.biojavax.ontology.ComparableTerm",
            "OBSOLETE"
            );
    public static final ChangeType DESCRIPTION = new ChangeType(
            "This term's description has changed",
            "org.biojavax.ontology.ComparableTerm",
            "DESCRIPTION"
            );
    public static final ChangeType RANKEDCROSSREF = new ChangeType(
            "This term's ranked crossrefs have changed",
            "org.biojavax.ontology.ComparableTerm",
            "RANKEDCROSSREF"
            );
    
    /**
     * Returns the (optional) identifier associated with this term.
     * @return the string identifier.
     */
    public String getIdentifier();
    
    /**
     * Sets the (optional) identifier associated with this term.
     * @param identifier the identifier to give the term. Null will unset it.
     * @throws ChangeVetoException if the identifier is unacceptable.
     */
    public void setIdentifier(String identifier) throws ChangeVetoException;
    
    /**
     * Checks to see if this term is obsolete. As the column in the database
     * is nullable, this value is a Boolean object instead of a boolean simple
     * type. Hence it may also be null.
     * @return true if it is, false if not.
     */
    public Boolean getObsolete();
    
    /**
     * Sets the obsolescence flag associated with this term.
     * @param obsolete true if it is obsolete, false if not. Nullable.
     * @throws ChangeVetoException if the change is unacceptable.
     */
    public void setObsolete(Boolean obsolete) throws ChangeVetoException;
    
    /**
     * Sets the description associated with this term.
     * @param description the description to give the term. Nullable.
     * @throws ChangeVetoException if the description is unacceptable.
     */
    public void setDescription(String description) throws ChangeVetoException;
    
}