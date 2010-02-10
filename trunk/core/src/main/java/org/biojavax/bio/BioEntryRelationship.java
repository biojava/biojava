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
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;
import org.biojavax.ontology.ComparableTerm;


/**
 * Represents the relation between two bioentries. The bioentry_relationship in
 * BioSQL is what this represents.
 * @author Mark Schreiber
 * @author Richard Holland
 * @see BioEntry
 * @since 1.5
 */
public interface BioEntryRelationship extends Comparable,Changeable {
    
    public static final ChangeType RANK = new ChangeType(
            "This bioentry relationship's rank has changed",
            "org.biojavax.bio.BioEntryRelationship",
            "RANK"
            );
    
    /**
     * Sets the rank of this relationship. The rank may be null in
     * the database, hence the use of an Integer object here and not
     * an int primitive.
     * @param rank Value of property rank.
     * @throws ChangeVetoException if the rank rankles.
     */
    public void setRank(Integer rank) throws ChangeVetoException;
    
    /**
     * Returns the rank of this relationship. The rank may be null in
     * the database, hence the use of an Integer object here and not
     * an int primitive.
     * @return Value of property rank.
     */
    public Integer getRank();
        
    /**
     * Returns the object of this relationship (ie. the BioEntry which
     * this relationship starts from). This is an immutable
     * property set by the constructor of an instantiating class.
     * @return Value of property object.
     */
    public BioEntry getObject();
    
    /**
     * Returns the subject of this relationship (ie. the BioEntry which
     * this relationship targets). This is an immutable
     * property set by the constructor of an instantiating class.
     * @return Value of property subject.
     */
    public BioEntry getSubject();
    
    /**
     * Returns the term describing the relationship. This is an immutable
     * property set by the constructor of an instantiating class.
     * @return Value of property term.
     */
    public ComparableTerm getTerm();
}

