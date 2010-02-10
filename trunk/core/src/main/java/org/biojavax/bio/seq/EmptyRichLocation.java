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
import java.util.Collections;
import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Unchangeable;
import org.biojavax.CrossRef;
import org.biojavax.CrossReferenceResolver;
import org.biojavax.RichAnnotation;
import org.biojavax.ontology.ComparableTerm;

/**
 * An Empty implementation of RichLocation. This class is intended to 
 * act as a place holder for events like the intersection of two locations
 * that do not overlap so that null need not be returned.
 * @author Richard Holland
 * @author Mark Schreiber
 * @author George Waldon
 * @since 1.5
 */
public class EmptyRichLocation extends Unchangeable implements RichLocation {
           
    /**
     * {@inheritDoc} 
     * ALWAYS RETURNS NULL
     */
    public RichFeature getFeature() { return null; }
        
    /**
     * {@inheritDoc} 
     * DOES NOTHING
     */
    public void sort() {}
    
    /**
     * {@inheritDoc} 
     * NOT IMPLEMENTED
     */
    public void setFeature(RichFeature feature) throws ChangeVetoException {
        throw new ChangeVetoException("Cannot set a feature for the empty location");
    }
    
    /**
     * {@inheritDoc} 
     * ALWAYS RETURNS NULL
     */
    public CrossRef getCrossRef() { return null; }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS THE EMPTY ANNOTATION
     */
    public Annotation getAnnotation() { return getRichAnnotation(); }

    /**
     * {@inheritDoc}
     * ALWAYS RETURNS THE EMPTY ANNOTATION
     */
    public RichAnnotation getRichAnnotation() { return RichAnnotation.EMPTY_ANNOTATION; }

    /**
     * {@inheritDoc} 
     * ALWAYS RETURNS THE EMPTY ANNOTATION NOTE SET
     */
    public Set getNoteSet() { return RichAnnotation.EMPTY_ANNOTATION.getNoteSet(); }
    
    /**
     * {@inheritDoc} 
     * NOT IMPLEMENTED
     */
    public void setNoteSet(Set notes) throws ChangeVetoException {
        throw new ChangeVetoException("Cannot annotate the empty location");
    }
    
    /**
     * {@inheritDoc} 
     * ALWAYS RETURNS NULL
     */
    public ComparableTerm getTerm() { return null; }
    
    /**
     * {@inheritDoc} 
     * NOT IMPLEMENTED
     */
    public void setTerm(ComparableTerm term) throws ChangeVetoException {
        throw new ChangeVetoException("Cannot give a term to the empty location");
    }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS ZERO
     */
    public int getCircularLength() { return 0; }
    
    /**
     * {@inheritDoc} 
     * NOT IMPLEMENTED
     */
    public void setCircularLength(int sourceSeqLength) throws ChangeVetoException {
        throw new ChangeVetoException("Cannot make empty locations circular");
    }
    
    /**
     * {@inheritDoc} 
     * ALWAYS RETURNS THE UNKNOWN STRAND
     */
    public Strand getStrand() { return Strand.UNKNOWN_STRAND; }
        
    /**
     * {@inheritDoc} 
     * ALWAYS RETURNS ZERO
     */
    public int getRank() { return 0; }
    
    /**
     * {@inheritDoc} 
     * NOT IMPLEMENTED
     */
    public void setRank(int rank) throws ChangeVetoException {
        throw new ChangeVetoException("Cannot give a rank to the empty location");
    }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS ZERO
     */
    public int getMax() { return 0; }
        
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS ZERO
     */
    public int getMin() { return 0; }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS THE EMPTY POSITION
     */ 
    public Position getMinPosition() { return Position.EMPTY_POSITION; }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS THE EMPTY POSITION
     */ 
    public Position getMaxPosition() { return Position.EMPTY_POSITION; }
    
    /**
     * {@inheritDoc} This method is ignored in the empty location because positions
     * are fixed an cannot be modified.
     */
    public void setPositionResolver(PositionResolver p) {} // ignore
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS THE EMPTY SET ITERATOR
     */
    public Iterator blockIterator() { return Collections.EMPTY_SET.iterator(); }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS TRUE
     */
    public boolean isContiguous() { return true; }
        
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS FALSE
     */
    public boolean contains(int p) { return false; }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS NULL
     */
    public Location getDecorator(Class decoratorClass) { return null; }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS PASSED LOCATION
     */
    public Location newInstance(Location loc) { return loc; }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS SELF
     */
    public Location translate(int dist) { return this; }  
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS FALSE
     */
    public boolean contains(Location l) { return false; }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS FALSE
     */
    public boolean overlaps(Location l) { return false; }
    
    /**
     * {@inheritDoc} 
     * ALWAYS RETURNS PASSED LOCATION
     */
    public Location union(Location l) {
        if (l==null) throw new IllegalArgumentException("Location cannot be null");
        if (!(l instanceof RichLocation)) l = RichLocation.Tools.enrich(l);
        return l;
    }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS SELF
     */
    public Location intersection(Location l) {
        if (l==null) throw new IllegalArgumentException("Location cannot be null");
        return this;
    }
       
    
    /**
     * {@inheritDoc} This method is ignored in the empty location because 
     * there is nothing to resolve.
     */
    public void setCrossRefResolver(CrossReferenceResolver r) {}
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS THE EMPTY SYMBOL LIST
     */
    public SymbolList symbols(SymbolList seq) {
        if (seq==null) throw new IllegalArgumentException("Sequence cannot be null");
        return SymbolList.EMPTY_LIST;
    }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS 17
     */
    public int hashCode() { return 17; }
    
    /**
     * {@inheritDoc}
     * Empty Rich Locations only match other Empty Rich Locations
     */
    public boolean equals(Object o) {
        if (o instanceof EmptyRichLocation) return true;
        return false;
    }
    
    /**
     * {@inheritDoc}
     * Empty Rich Locations return 0 when compared to other Empty ones,
     * or -1 otherwise.
     */
    public int compareTo(Object o) {
        if (o instanceof EmptyRichLocation) return 0;
        else return -1;
    }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS "empty"
     */
    public String toString() {
        return "empty";
    }
}

