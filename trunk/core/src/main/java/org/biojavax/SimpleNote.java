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
import org.biojavax.ontology.ComparableTerm;

/**
 * Simple implementation of Note.
 * @author Richard Holland
 * @author George Waldon - limited firing
 * @since 1.5
 */
public class SimpleNote extends AbstractChangeable implements Note {
    
    private ComparableTerm term;
    private String value;
    private int rank;
    
    /**
     * Creates a new instance of SimpleNote with a given term, value and rank.
     * @param term the term of the note. Cannot be null.
     * @param value the (optional) value to give it.
     * @param rank the rank to give it.
     */
    public SimpleNote(ComparableTerm term, String value, int rank) {
        if (term==null) throw new IllegalArgumentException("Term cannot be null");
        this.term = term;
        this.value = value;
        this.rank = rank;
    }
    
    // Hibernate requirement - not for public use.
    protected SimpleNote() {}
    
    /**
     * {@inheritDoc}
     */
    public ComparableTerm getTerm() { return this.term; }
    
    /**
     * {@inheritDoc}
     */
    public void setTerm(ComparableTerm term) throws ChangeVetoException {
        if (term==null) throw new IllegalArgumentException("Term cannot be null");
        if(term.equals(this.term))
            return;
        if(!this.hasListeners(Note.TERM)) {
            this.term = term;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    Note.TERM,
                    term,
                    this.term
                    );
            ChangeSupport cs = this.getChangeSupport(Note.TERM);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.term = term;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public String getValue() { return this.value; }
    
    /**
     * {@inheritDoc}
     */
    public void setValue(String value) throws ChangeVetoException {
        if(this.value!=null && this.value.equals(value)) return;
        else if(this.value==null && value==null) return;
        
        if(!this.hasListeners(Note.VALUE)) {
            this.value = value;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    Note.VALUE,
                    value,
                    this.value
                    );
            ChangeSupport cs = this.getChangeSupport(Note.VALUE);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.value = value;
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
    public void setRank(int rank) throws ChangeVetoException {
        if(this.rank==rank) 
            return;
        if(!this.hasListeners(Note.RANK)) {
            this.rank = rank;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    Note.RANK,
                    new Integer(rank),
                    new Integer(this.rank)
                    );
            ChangeSupport cs = this.getChangeSupport(Note.RANK);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.rank = rank;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     * Notes are compared first by rank, then by the term.
     */
    public int compareTo(Object o) {
        if (o==this) return 0;
        // Hibernate comparison - we haven't been populated yet
        if (this.term==null) return -1;
        // Normal comparison
        Note them = (Note)o;
        if (this.rank!=them.getRank()) return this.rank-them.getRank();
        else return this.term.compareTo(them.getTerm());
    }
    
    /**
     * {@inheritDoc}
     * Notes are equal if they have the same rank and term.
     */
    public boolean equals(Object o) {
        if (o==this) return true;
        if (!(o instanceof Note)) return false;
        // Hibernate comparison - we haven't been populated yet
        if (this.term==null) return false;
        // Normal comparison
        Note them = (Note)o;
        return this.term.equals(them.getTerm()) && this.rank==them.getRank();
    }
    
    /**
     * {@inheritDoc}
     */
    public int hashCode() {
        int hash = 17;
        // Hibernate comparison - we haven't been populated yet
        if (this.term==null) return hash;
        // Normal comparison
        hash = 31*hash + this.term.hashCode();
        hash = 31*hash + this.rank;
        return hash;
    }
    
    /**
     * {@inheritDoc}
     * Form: "(#rank) term: value"
     */
    public String toString() {
        return "(#"+this.rank+") "+this.term+": "+this.value;
    }
}
