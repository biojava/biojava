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

import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;

/**
 * Simple implementation of RankedCrossRef.
 * @author Richard Holland
 * @author gwaldon
 * @since 1.5
 */
public class SimpleRankedCrossRef extends AbstractChangeable implements RankedCrossRef {
    
    private CrossRef crossref;
    private int rank;
    
    /**
     * Constructs a new crossref with a rank.
     * @param crossref the crossref to rank. Must not be null.
     * @param rank the rank to give it.
     */
    public SimpleRankedCrossRef(CrossRef crossref, int rank) {
        if (crossref==null) throw new IllegalArgumentException("Cross reference cannot be null");
        this.crossref = crossref;
        this.rank = rank;
    }
    
    // Hibernate requirement - not for public use.
    protected SimpleRankedCrossRef() {}
    
    // Hibernate requirement - not for public use.
    void setCrossRef(CrossRef crossref) { this.crossref = crossref; }
    
    /**
     * {@inheritDoc}
     */
    public CrossRef getCrossRef() { return this.crossref; }
    
    /**
     * {@inheritDoc}
     */
    public void setRank(int rank)  throws ChangeVetoException {
        if(rank==this.rank)
            return;
        if(!this.hasListeners(RankedCrossRef.RANK)) {
            this.rank = rank;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    RankedCrossRef.RANK,
                    new Integer(rank),
                    new Integer(this.rank)
                    );
            ChangeSupport cs = this.getChangeSupport(RankedCrossRef.RANK);
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
     * Ranked cross references are the same if they have the same rank and
     * refer to the same cross reference (cross references are equal).
     */
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj==null || !(obj instanceof RankedCrossRef)) return false;
        // Hibernate comparison - we haven't been populated yet
        if (this.crossref==null) return false;
        // Normal comparison
        RankedCrossRef them = (RankedCrossRef)obj;
        return (this.rank==them.getRank() &&
                this.crossref.equals(them.getCrossRef()));
    }
    
    /**
     * {@inheritDoc}
     * Ranked cross references are sorted first by rank, then by cross reference.
     */
    public int compareTo(Object o) {
        if (o==this) return 0;
        // Hibernate comparison - we haven't been populated yet
        if (this.crossref==null) return -1;
        // Normal comparison
        RankedCrossRef them = (RankedCrossRef)o;
        if (this.rank!=them.getRank()) return this.rank - them.getRank();
        return this.crossref.compareTo(them.getCrossRef());
    }
    
    /**
     * {@inheritDoc}
     */
    public int hashCode() {
        int code = 17;
        // Hibernate comparison - we haven't been populated yet
        if (this.crossref==null) return code;
        // Normal comparison
        code = 37*code + this.crossref.hashCode();
        code = 37*code + this.rank;
        return code;
    }
    
    /**
     * {@inheritDoc}
     * Form: "(#rank) crossref"
     */
    public String toString() {
        return "(#"+this.rank+") "+this.crossref;
    }
}
