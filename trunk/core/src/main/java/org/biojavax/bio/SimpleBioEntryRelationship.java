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

/*
 * SimpleBioEntryRelationship.java
 *
 * Created on June 16, 2005, 2:07 PM
 */

package org.biojavax.bio;

import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.ontology.ComparableTerm;

/**
 * Represents a relationship between two bioentries that is described by a term 
 * and given a rank.
 * @author Richard Holland
 * @author Mark Schreiber
 * @author George Waldon
 * @since 1.5
 */
public class SimpleBioEntryRelationship extends AbstractChangeable implements BioEntryRelationship {
    
    private BioEntry object;
    private BioEntry subject;
    private ComparableTerm term;
    private Integer rank;
    
    /**
     * Creates a new instance of SimpleBioEntryRelationship. None of the parameters
     * may be null, and all are immutable except the rank.
     * @param rank The rank of the relationship.
     * @param subject The subject bioentry.
     * @param term The relationship term.
     */
    public SimpleBioEntryRelationship(BioEntry object, BioEntry subject, ComparableTerm term, Integer rank) {
        if (object==null) throw new IllegalArgumentException("Object cannot be null");
        if (subject==null) throw new IllegalArgumentException("Subject cannot be null");
        if (term==null) throw new IllegalArgumentException("Term cannot be null");
        this.object = object;
        this.subject = subject;
        this.term = term;
        this.rank = rank;
    }
    
    // Hibernate requirement - not for public use.
    protected SimpleBioEntryRelationship() {}
    
    /**
     * {@inheritDoc}
     */
    public void setRank(Integer rank) throws ChangeVetoException {
        if(!this.hasListeners(BioEntryRelationship.RANK)) {
            this.rank = rank;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    BioEntryRelationship.RANK,
                    rank,
                    this.rank
                    );
            ChangeSupport cs = this.getChangeSupport(BioEntryRelationship.RANK);
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
    public Integer getRank() { return this.rank; }
        
    /**
     * {@inheritDoc}
     */
    public BioEntry getObject() { return this.object; }
    
    // Hibernate requirement - not for public use.
    void setObject(BioEntry object) { this.object = object; }
        
    /**
     * {@inheritDoc}
     */
    public BioEntry getSubject() { return this.subject; }
    
    // Hibernate requirement - not for public use.
    void setSubject(BioEntry subject) { this.subject = subject; }
    
    /**
     * {@inheritDoc}
     */
    public ComparableTerm getTerm() { return this.term; }
    
    // Hibernate requirement - not for public use.
    void setTerm(ComparableTerm term) { this.term = term; }
    
    /**
     * {@inheritDoc}
     * A relationship is compared first by rank, then object, subject, and term. If
     * ranks are null, they are treated as zero for comparison's sake.
     */
    public int compareTo(Object o) {
        if (o==this) return 0;
        // Hibernate comparison - we haven't been populated yet
        if (this.object==null) return -1;
        // Normal comparison
        BioEntryRelationship them = (BioEntryRelationship)o;
        int ourRank =   (this.getRank() == null ? 0 : this.rank.intValue());
        int theirRank = (them.getRank() == null ? 0 : them.getRank().intValue());
        if (ourRank!=theirRank) return ourRank-theirRank;
        if (!this.object.equals(them.getObject())) return this.object.compareTo(them.getObject());
        if (!this.subject.equals(them.getSubject())) return this.subject.compareTo(them.getSubject());
        return this.term.compareTo(them.getTerm());
    }
    
    /**
     * {@inheritDoc} 
     * Relationships are equal if they share the same rank, object, subject and term. If
     * ranks are null, they are treated as zero for consistency with comparison.
     */
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj==null || !(obj instanceof BioEntryRelationship)) return false;
        // Hibernate comparison - we haven't been populated yet
        if (this.object==null) return false;
        // Normal comparison
        BioEntryRelationship them = (BioEntryRelationship)obj;
        int ourRank =   (this.getRank() == null ? 0 : this.rank.intValue());
        int theirRank = (them.getRank() == null ? 0 : them.getRank().intValue());
        return (this.subject.equals(them.getSubject()) &&
                    this.object.equals(them.getObject()) &&
                    this.term.equals(them.getTerm()) &&
                    ourRank==theirRank);
    }
    
    /**
     * {@inheritDoc}
     */
    public int hashCode() {
        int code = 17;
        // Hibernate comparison - we haven't been populated yet
        if (this.subject==null) return code;
        // Normal comparison
        code = code*37 + (this.getRank() == null ? 0 : this.rank.hashCode());
        code = code*37 + this.object.hashCode();
        code = code*37 + this.subject.hashCode();
        code = code*37 + this.term.hashCode();
        return code;
    }
    
    /**
     * {@inheritDoc}
     * Form is "(#rank) term(object,subject)"
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

