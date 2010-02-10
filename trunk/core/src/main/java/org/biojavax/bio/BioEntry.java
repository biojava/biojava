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

package org.biojavax.bio;

import java.util.Set;

import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;
import org.biojavax.Comment;
import org.biojavax.Namespace;
import org.biojavax.RankedCrossRefable;
import org.biojavax.RankedDocRef;
import org.biojavax.RichAnnotatable;
import org.biojavax.bio.taxa.NCBITaxon;

/**
 * This class relates to the bioentry table in BioSQL. It holds everything you need
 * to define a non-sequence bearing bioentry.
 * @author Mark Schreiber
 * @author Richard Holland
 * @see BioEntryRelationship
 * @since 1.5
 */
public interface BioEntry extends RichAnnotatable,RankedCrossRefable,Comparable,Changeable {
    
    public static final ChangeType IDENTIFIER = new ChangeType(
            "This bioentry's identifier has changed",
            "org.biojavax.bio.BioEntry",
            "IDENTIFIER"
            );
    public static final ChangeType DESCRIPTION = new ChangeType(
            "This bioentry's description has changed",
            "org.biojavax.bio.BioEntry",
            "DESCRIPTION"
            );
    public static final ChangeType DIVISION = new ChangeType(
            "This bioentry's division has changed",
            "org.biojavax.bio.BioEntry",
            "DIVISION"
            );
    public static final ChangeType TAXON = new ChangeType(
            "This bioentry's taxon has changed",
            "org.biojavax.bio.BioEntry",
            "TAXON"
            );
    public static final ChangeType SEQVERSION = new ChangeType(
            "This bioentry's sequence version has changed",
            "org.biojavax.bio.BioEntry",
            "SEQVERSION"
            );
    public static final ChangeType RANKEDCROSSREF = new ChangeType(
            "This bioentry's ranked crossrefs changed",
            "org.biojavax.bio.BioEntry",
            "RANKEDCROSSREF"
            );
    public static final ChangeType RANKEDDOCREF = new ChangeType(
            "This bioentry's ranked docrefs changed",
            "org.biojavax.bio.BioEntry",
            "RANKEDDOCREF"
            );
    public static final ChangeType COMMENT = new ChangeType(
            "This bioentry's comments changed",
            "org.biojavax.bio.BioEntry",
            "COMMENT"
            );
    public static final ChangeType RELATIONS = new ChangeType(
            "This bioentry's relations have changed",
            "org.biojavax.bio.BioEntry",
            "RELATIONS"
            );
    
    /**
     * Returns the namespace of this bioentry. The namespace is supposed
     * to be an immutable property set by the constructor.
     * @return the namespace of this bioentry.
     */
    public Namespace getNamespace();
    
    /**
     * Returns the name of this bioentry. The name is supposed
     * to be an immutable property set by the constructor.
     * @return Value of property name.
     */
    public String getName();
    
    /**
     * Returns the accession of this bioentry. The accession is supposed
     * to be an immutable property set by the constructor.
     * @return Value of property accession.
     */
    public String getAccession();
    
    /**
     * Returns the identifier of this bioentry.
     * @return Value of property identifier.
     */
    public String getIdentifier();
    
    /**
     * Sets the identifier of this bioentry. Null is allowable.
     * @param identifier New value of property identifier.
     * @throws ChangeVetoException in case of objections.
     */
    public void setIdentifier(String identifier) throws ChangeVetoException;
    
    /**
     * Returns the division of this bioentry. Division relates to a division
     * of the parent namespace.
     * @return Value of property division.
     */
    public String getDivision();
    
    /**
     * Sets the division of this bioentry. Null is allowable.
     * @param division New value of property division.
     * @throws ChangeVetoException in case of objections.
     */
    public void setDivision(String division) throws ChangeVetoException;
    
    /**
     * Returns a description of this sequence.
     * @return Value of property description.
     */
    public String getDescription();
    
    /**
     * Sets the description for this bioentry.
     * @param description New value of property description.
     * @throws ChangeVetoException in case of objections.
     */
    public void setDescription(String description) throws ChangeVetoException;
    
    /**
     * Gets the version of this bioentry. Bioentries with no versions return 0.
     * The version is supposed to be immutable and set only by the constructor.
     * @return Value of property version.
     */
    public int getVersion();
    
    /**
     * Gets the taxon associated with this bioentry. It may be null.
     * @return Value of property taxon.
     */
    public NCBITaxon getTaxon();
    
    /**
     * Sets the taxon for this bioentry. It may be null, in which case the
     * taxon is unset.
     * @param taxon New value of property taxon.
     * @throws ChangeVetoException in case of objections.
     */
    public void setTaxon(NCBITaxon taxon) throws ChangeVetoException;
    
    /**
     * Returns a set of all bioentrydocrefs associated with this bioentry. This
     * set is not mutable. If no docrefs are associated, you will get back an
     * empty set.
     * @return a set of RankedDocRef objects.
     * @see RankedDocRef
     */
    public Set getRankedDocRefs();
    
    /**
     * Returns a set of all comments associated with this bioentry. This
     * set is not mutable. If no comments are associated, you will get back an
     * empty set.
     * @return a set of Comment objects.
     * @see Comment
     */
    public Set getComments();
    
    /**
     * Returns a set of all relationships associated with this bioentry. This
     * set is not mutable. If no relationships are associated, you will get back an
     * empty set.
     * @return a set of BioEntryRelationship objects.
     * @see BioEntryRelationship
     */
    public Set getRelationships();
    
    /**
     * Adds a ranked docref instance to this bioentry. Must not be null.
     * @param docref the item to add.
     * @throws ChangeVetoException if it doesn't want to add it.
     */
    public void addRankedDocRef(RankedDocRef docref) throws ChangeVetoException;
    
    /**
     * Removes a ranked docref instance from this bioentry. If it was not found,
     * nothing happens.
     * @param docref the item to remove.
     * @throws ChangeVetoException if it doesn't want to remove it.
     */
    public void removeRankedDocRef(RankedDocRef docref) throws ChangeVetoException;
    
    /**
     * Adds a comment instance to this bioentry. Must not be null.
     * @param comment the item to add.
     * @throws ChangeVetoException if it doesn't want to add it.
     */
    public void addComment(Comment comment) throws ChangeVetoException;
    
    /**
     * Removes a comment instance from this bioentry. If it wasn't present, it
     * nothing will happen.
     * @param comment the item to remove.
     * @throws ChangeVetoException if it doesn't want to remove it.
     */
    public void removeComment(Comment comment) throws ChangeVetoException;
    
    /**
     * Adds a relation instance to this bioentry. It must not be null.
     * @param relation the item to add.
     * @throws ChangeVetoException if it doesn't want to add it.
     */
    public void addRelationship(BioEntryRelationship relation) throws ChangeVetoException;
    
    /**
     * Removes a relation instance from this bioentry. If it wasn't present,
     * nothing will happen.
     * @param relation the item to remove.
     * @throws ChangeVetoException if it doesn't want to remove it.
     */
    public void removeRelationship(BioEntryRelationship relation) throws ChangeVetoException;
}



