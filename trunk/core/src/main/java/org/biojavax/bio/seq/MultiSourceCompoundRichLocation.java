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
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.symbol.Location;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.CrossRef;
import org.biojavax.ontology.ComparableTerm;

/**
 * An implementation of RichLocation which possibly covers multiple locations, 
 * on different strands, different circular lengths, or different sequences.
 * @author Richard Holland
 * @author Mark Schreiber
 * @since 1.5
 */
public class MultiSourceCompoundRichLocation extends CompoundRichLocation implements RichLocation {
    
    /**
     * Constructs a MultiSourceCompoundRichLocation from the given set of members, with
     * the default term of "join". Note that you really shouldn't use this if
     * you are unsure if your members set contains overlapping members. Use
     * RichLocation.Tools.construct() instead. The members collection
     * must only contain Location instances. Any that are not RichLocations will
     * be converted using RichLocation.Tools.enrich().
     * @param members the members to put into the compound location.
     * @see RichLocation.Tools
     */
    public MultiSourceCompoundRichLocation(Collection members) { this(getJoinTerm(), members); }
    
    /**
     * Constructs a MultiSourceCompoundRichLocation from the given set of members.
     * Note that you really shouldn't use this if
     * you are unsure if your members set contains overlapping members. Use
     * RichLocation.Tools.construct(members) instead. The members collection
     * must only contain Location instances. Any that are not RichLocations will
     * be converted using RichLocation.Tools.enrich().
     * @param term the term to use when describing the group of members.
     * @param members the members to put into the compound location.
     * @see RichLocation.Tools
     */
    public MultiSourceCompoundRichLocation(ComparableTerm term, Collection members) {
        if (term==null) throw new IllegalArgumentException("Term cannot be null");
        if (members==null || members.size()<2) throw new IllegalArgumentException("Must have at least two members");        
        this.term = term;
        this.members = new ArrayList();   
        for (Iterator i = members.iterator(); i.hasNext(); ) {
            // Convert each member into a RichLocation
            Object o = i.next();
            if (!(o instanceof RichLocation)) o = RichLocation.Tools.enrich((Location)o);
            // Convert
            RichLocation rl = (RichLocation)o;
            // Add in member
            this.members.add(rl);
            // Update our size
            this.size += Math.max(rl.getMin(),rl.getMax())-Math.min(rl.getMin(),rl.getMax());
        }
    }
                   
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS NULL
     */
    public CrossRef getCrossRef() { return null; }
            
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS ZERO
     */
    public int getCircularLength() { return 0; }
    
    /**
     * {@inheritDoc}
     * NOT IMPLEMENTED
     * @throws ChangeVetoException ALWAYS
     */
    public void setCircularLength(int sourceSeqLength) throws ChangeVetoException {
        if (sourceSeqLength>0) throw new ChangeVetoException("MultiSourceCompoundRichLocations cannot be circular");
    }
    
    /**
     * {@inheritDoc}
     */
    public Strand getStrand() { return Strand.UNKNOWN_STRAND; }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS ONE
     */
    public int getMin() { return 1; }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS COMBINED LENGTH OF MEMBERS
     */
    public int getMax() { return this.size; }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS A POINT POSITION AT POINT 1
     */
    public Position getMinPosition() { return new SimplePosition(false,false,1); }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS A POINT POSITION AT POINT EQUIVALENT TO COMBINED LENGTH OF MEMBERS
     */
    public Position getMaxPosition() { return new SimplePosition(false,false,this.size); }
    
    /**
     * {@inheritDoc}
     * Recursively translates all members of this location.
     */
    public Location translate(int dist) {
        if (this.members.isEmpty()) return this;
        List newmembers = new ArrayList();
        for (Iterator i = this.members.iterator(); i.hasNext(); ) {
            RichLocation rl = (RichLocation)i.next();
            newmembers.add(rl.translate(dist));
        }
        return new MultiSourceCompoundRichLocation(this.getTerm(),newmembers);
    }    
}

