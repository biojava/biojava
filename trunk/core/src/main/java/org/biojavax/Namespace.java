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

import java.net.URI;

import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * The namespace of an entry in a database schema. Relates directly to the
 * BioDatabase table in BioSQL. All BioEntry objects belong to namespaces.
 * @author Mark Schreiber
 * @author Richard Holland
 * @see org.biojavax.bio.BioEntry
 * @since 1.5
 */
public interface Namespace extends Comparable,Changeable {
    
    public static final ChangeType NAME = new ChangeType(
            "This namespace's name has changed",
            "org.biojavax.Namespace",
            "NAME"
            );
    public static final ChangeType AUTHORITY = new ChangeType(
            "This namespace's authority has changed",
            "org.biojavax.Namespace",
            "AUTHORITY"
            );
    public static final ChangeType DESCRIPTION = new ChangeType(
            "This namespace's description has changed",
            "org.biojavax.Namespace",
            "DESCRIPTION"
            );
    public static final ChangeType ACRONYM = new ChangeType(
            "This namespace's acronym has changed",
            "org.biojavax.Namespace",
            "ACRONYM"
            );
    public static final ChangeType URI = new ChangeType(
            "This namespace's URI has changed",
            "org.biojavax.Namespace",
            "URI"
            );
    
    /**
     * The name of the namespace is immutable and must be set by the constructor
     * of the instantiating class. The name should also be unique. This method
     * will return the name.
     * @return The name of the namespace.
     */
    public String getName();
    
    /**
     * This method will return the authority that governs the namespace.
     * @return the name of the namespace authority.
     */
    public String getAuthority();
    
    /**
     * This method sets the authority that governs the namespace. Null will 
     * unset it. 
     * @param authority the name of the namespace authority.
     * @throws ChangeVetoException in case of objections.
     */
    public void setAuthority(String authority) throws ChangeVetoException;
    
    /**
     * Returns a description of this namespace.
     * @return the description of the namespace.
     */
    public String getDescription();
    
    /**
     * This method sets a description for the namespace. Null will unset it. 
     * @param description the description of the namespace.
     * @throws ChangeVetoException in case of objections.
     */
    public void setDescription(String description) throws ChangeVetoException;
    
    /**
     * If the namespace has an acronym, this will return it.
     * @return the acronym for the namespace.
     */
    public String getAcronym();
    
    /**
     * Sets an optional acronym for the namespace. Null will unset it. Note that
     * in BioSQL 1.0 Acronym is only part of the Oracle schema therefore it will
     * only be persisted in that schema.
     * @param acronym the acronym for the namespace.
     * @throws ChangeVetoException in case of objections.
     */
    public void setAcronym(String acronym) throws ChangeVetoException;
    
    /**
     * If the namespace has a URI, this will return it.
     * @return the URI of the authority.
     */
    public URI getURI();
    
    /**
     * Sets an optional URI for the namespace. Null will unset it. Note that in
     * BioSQL 1.0 URI is not persisted into the database unless the 
     * extended Oracle schema is used.
     * @param URI the URI of the authority.
     * @throws ChangeVetoException in case of objections.
     */
    public void setURI(URI URI) throws ChangeVetoException;
    
}

