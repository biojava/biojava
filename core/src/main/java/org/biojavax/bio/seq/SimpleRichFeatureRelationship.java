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

import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.RichObjectFactory;
import org.biojavax.ontology.ComparableTerm;

/**
 * Represents a relationship between two features that is described by a term.
 * @author Richard Holland
 * @author Mark Schreiber
 * @since 1.5
 *
 */
public class SimpleRichFeatureRelationship extends AbstractChangeable implements RichFeatureRelationship {

    private RichFeature object;
    private RichFeature subject;
    private ComparableTerm term;
    private int rank;
    
    /**
     * Gets the default CONTAINS term used for defining the relationship between features.
     * @return the default CONTAINS term.
     */
    public static ComparableTerm getContainsTerm() {
        return RichObjectFactory.getDefaultOntology().getOrCreateTerm("contains");
    }
    
    /**
     * Creates a new instance of SimpleRichFeatureRelationship.
     * @param subject The subject RichFeature.
     * @param term The relationship term.
     * @param rank the rank of the relationship.
     */    
    public SimpleRichFeatureRelationship(RichFeature object, RichFeature subject, ComparableTerm term, int rank) {
        if (object==null) throw new IllegalArgumentException("Object cannot be null");
        if (subject==null) throw new IllegalArgumentException("Subject cannot be null");
        if (term==null) throw new IllegalArgumentException("Term cannot be null");
        this.object = object;
        this.subject = subject;
        this.term = term;
        this.rank = rank;
    }
    
    // Hibernate requirement - not for public use.
    protected SimpleRichFeatureRelationship() {}
    
    /**
     * {@inheritDoc}
     */
    public void setRank(int rank) throws ChangeVetoException {
        if(!this.hasListeners(RichFeatureRelationship.RANK)) {
            this.rank = rank;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    RichFeatureRelationship.RANK,
                    new Integer(rank),
                    new Integer(this.rank)
                    );
            ChangeSupport cs = this.getChangeSupport(RichFeatureRelationship.RANK);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.rank = rank;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public int getRank() { return this.rank; }
        
    /**
     * {@inheritDoc}
     */
    public RichFeature getObject() { return this.object; }
    
    // Hibernate requirement - not for public use.
    void setObject(RichFeature object) { this.object = object; }
        
    /**
     * {@inheritDoc}
     */
    public RichFeature getSubject() { return this.subject; }
    
    // Hibernate requirement - not for public use.
    void setSubject(RichFeature subject) { this.subject = subject; }
    
    /**
     * {@inheritDoc}
     */
    public ComparableTerm getTerm() { return this.term; }
    
    // Hibernate requirement - not for public use.
    void setTerm(ComparableTerm term) { this.term = term; }
    
    /**
     * {@inheritDoc}
     * Relations are compared first by rank, then object, subject, then finally term.
     */
    public int compareTo(Object o) {
        if (o==this) return 0;
        // Hibernate comparison - we haven't been populated yet
        if (this.object==null) return -1;
        // Normal comparison
        RichFeatureRelationship them = (RichFeatureRelationship)o;
        if (this.rank!=them.getRank()) return this.rank-them.getRank();
        if (!this.object.equals(them.getObject())) return this.object.compareTo(them.getObject());
        if (!this.subject.equals(them.getSubject())) return this.subject.compareTo(them.getSubject());
        else return this.getTerm().compareTo(them.getTerm());
    }
    
    /**
     * {@inheritDoc}
     * Relations are equal if their objects, subjects and terms are equal.
     */
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj==null || !(obj instanceof RichFeatureRelationship)) return false;
        // Hibernate comparison - we haven't been populated yet
        if (this.object==null) return false;
        // Normal comparison
        RichFeatureRelationship them = (RichFeatureRelationship)obj;
        return (this.object.equals(them.getObject()) &&
                this.subject.equals(them.getSubject()) &&
                this.term.equals(them.getTerm()));
    }
    
    /**
     * {@inheritDoc}
     */
    public int hashCode() {
        int code = 17;
        // Hibernate comparison - we haven't been populated yet
        if (this.subject==null) return code;
        // Normal comparison
        code = code*37 + this.object.hashCode();
        code = code*37 + this.subject.hashCode();
        code = code*37 + this.term.hashCode();
        return code;
    }
    
    /**
     * {@inheritDoc}
     * Form: "(#rank) term(object,subject)"
     */
    public String toString() {
        return "(#"+this.rank+") "+this.getTerm()+"("+this.getObject()+","+this.getSubject()+")";
    }
    
    // Hibernate requirement - not for public use.
    private Integer id;
    
    /**
     * Gets the Hibernate ID. Should be used with caution.
     * @return the Hibernate ID, if using Hibernate.
     */
    public Integer getId() { return this.id; }
    
    /**
     * Sets the Hibernate ID. Should be used with caution.
     * @param id the Hibernate ID, if using Hibernate.
     */
    public void setId(Integer id) { this.id = id;}
}

