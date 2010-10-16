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
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;
import org.biojavax.ontology.ComparableTerm;

/**
 * Represents the relation between two features. The seqfeature_relationship 
 * in BioSQL is what this is based on.
 * @author Mark Schreiber
 * @author Richard Holland
 * @since 1.5
 */
public interface RichFeatureRelationship extends Comparable,Changeable {
        
    public static final ChangeType RANK = new ChangeType(
            "This feature relationship's rank has changed",
            "org.biojavax.bio.seq.RichFeatureRelationship",
            "RANK"
            );
    
    /**
     * Sets the rank of this relationship.
     * @param rank Value of property rank.
     * @throws ChangeVetoException if the rank is untasty.
     */
    public void setRank(int rank) throws ChangeVetoException;
    
    /**
     * Gets the rank of this relationship.
     * @return Value of property rank.
     */
    public int getRank();
        
    /**
     * Returns the object of this relationship (ie. the feature which
     * this relationship starts from). This is an immutable
     * property set by the constructor of an instantiating class.
     * @return Value of property object.
     */
    public RichFeature getObject();
    
    /**
     * Gets the feature that this relationship refers to. This is set
     * at constructor time and is immutable.
     * @return Value of property subject.
     */
    public RichFeature getSubject();
    
    /**
     * Gets the term that describes this relationship. This is set
     * at constructor time and is immutable.
     * @return Value of property term.
     */
    public ComparableTerm getTerm();
    
}
