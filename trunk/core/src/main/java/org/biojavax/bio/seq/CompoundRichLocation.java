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
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.biojava.bio.Annotation;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.CrossReferenceResolver;
import org.biojavax.RichAnnotation;
import org.biojavax.RichObjectFactory;
import org.biojavax.ontology.ComparableTerm;

/**
 * An implementation of RichLocation which covers multiple locations, 
 * but on the same strand of the same (optionally circular) sequence.
 * @author Richard Holland
 * @author Mark Schreiber
 * @since 1.5
 */
public class CompoundRichLocation extends SimpleRichLocation implements RichLocation {

    protected List members;
    protected int size = 0;

    /**
     * Getter for the "join" term
     * @return the "join" term
     */
    public static ComparableTerm getJoinTerm() {
        return RichObjectFactory.getDefaultOntology().getOrCreateTerm("join");
    }

    /**
     * Getter for the "order" term
     * @return the "order" term
     */
    public static ComparableTerm getOrderTerm() {
        return RichObjectFactory.getDefaultOntology().getOrCreateTerm("order");
    }

    /**
     * Constructs a CompoundRichLocation from the given set of members, with
     * the default term of "join". Note that you really shouldn't use this if
     * you are unsure if your members set contains overlapping members. Use
     * RichLocation.Tools.construct() instead. The members collection
     * must only contain Location instances. Any that are not RichLocations will
     * be converted using RichLocation.Tools.enrich(). All members must come from
     * the same strand of the same sequence with the same circular length.
     * @param members the members to put into the compound location.
     * @see RichLocation.Tools
     */
    public CompoundRichLocation(Collection members) {
        this(getJoinTerm(), members);
    }

    /**
     * Constructs a CompoundRichLocation from the given set of members.
     * Note that you really shouldn't use this if
     * you are unsure if your members set contains overlapping members. Use
     * RichLocation.Tools.construct(members) instead. The members collection
     * must only contain Location instances. Any that are not RichLocations will
     * be converted using RichLocation.Tools.enrich().
     * @param term the term to use when describing the group of members.
     * @param members the members to put into the compound location.
     * @see RichLocation.Tools
     */
    public CompoundRichLocation(ComparableTerm term, Collection members) {
        if (term == null) {
            throw new IllegalArgumentException("Term cannot be null");
        }
        if (members == null || members.size() < 2) {
            throw new IllegalArgumentException("Must have at least two members");
        }
        if (RichLocation.Tools.isMultiSource(members)) {
            throw new IllegalArgumentException("All members must be from the same source");
        }
        this.term = term;
        this.members = new ArrayList();
        for (Iterator i = members.iterator(); i.hasNext();) {
            // Convert each member into a RichLocation
            Object o = i.next();
            if (!(o instanceof RichLocation)) {
                o = RichLocation.Tools.enrich((Location) o);
            }
            // Convert
            RichLocation rl = (RichLocation) o;
            // Add in member
            this.members.add(rl);
            // Update our cross ref
            this.setCrossRef(rl.getCrossRef());
            // Update our circular length
            this.circularLength = rl.getCircularLength();
            // Update our strand
            this.setStrand(rl.getStrand());
            // Update our size and min/max
            this.size += rl.getMax() - rl.getMin();
            if (this.getMinPosition() == null) {
                this.setMinPosition(rl.getMinPosition());
            } else {
                this.setMinPosition(this.posmin(this.getMinPosition(), rl.getMinPosition()));
            }
            if (this.getMaxPosition() == null) {
                this.setMaxPosition(rl.getMaxPosition());
            } else {
                this.setMaxPosition(this.posmax(this.getMaxPosition(), rl.getMaxPosition()));
            }
        }
    }

    // for internal use only
    protected CompoundRichLocation() {
    }

    /**
     * {@inheritDoc} 
     */
    public void sort() {
        Collections.sort(this.members);
    }

    /**
     * {@inheritDoc} 
     * Passes the call on to each of its members in turn.
     */
    public void setFeature(RichFeature feature) throws ChangeVetoException {
        super.setFeature(feature);
        for (Iterator i = this.members.iterator(); i.hasNext();) {
            ((RichLocation) i.next()).setFeature(feature);
        }
    }

    /**
     * {@inheritDoc}
     * ALWAYS RETURNS THE EMPTY ANNOTATION
     */
    public Annotation getAnnotation() {
        return getRichAnnotation();
    }

    /**
     * {@inheritDoc}
     * ALWAYS RETURNS THE EMPTY ANNOTATION
     */
    public RichAnnotation getRichAnnotation() {
        return RichAnnotation.EMPTY_ANNOTATION;
    }
    
    /**
     * {@inheritDoc}
     * ALWAYS RETURNS THE EMPTY ANNOTATION NOTE SET
     */
    public Set getNoteSet() {
        return RichAnnotation.EMPTY_ANNOTATION.getNoteSet();
    }

    /**
     * {@inheritDoc}
     * NOT IMPLEMENTED
     * @throws ChangeVetoException ALWAYS
     */
    public void setNoteSet(Set notes) throws ChangeVetoException {
        throw new ChangeVetoException("Cannot annotate compound locations.");
    }

    /**
     * {@inheritDoc}
     * RECURSIVELY APPLIES CALL TO ALL MEMBERS
     */
    public void setCircularLength(int sourceSeqLength) throws ChangeVetoException {
        super.setCircularLength(sourceSeqLength);
        for (Iterator i = this.members.iterator(); i.hasNext();) {
            ((RichLocation) i.next()).setCircularLength(sourceSeqLength);
        }
    }

    /**
     * {@inheritDoc}
     */
    public Iterator blockIterator() {
        final List sortedMembers = new ArrayList(this.members);
        Collections.sort(sortedMembers);
        return sortedMembers.iterator();
    }

    /**
     * {@inheritDoc}
     * ALWAYS RETURNS FALSE
     */
    public boolean isContiguous() {
        return false;
    }

    /**
     * {@inheritDoc}
     * Recursively applies this call to all members.
     */
    public boolean contains(int p) {
        for (Iterator i = this.members.iterator(); i.hasNext();) {
            if (((RichLocation) i.next()).contains(p)) {
                return true;
            }
        }
        return false;
    }

    /**
     * {@inheritDoc}
     * ALWAYS RETURNS NULL
     */
    public Location getDecorator(Class decoratorClass) {
        return null;
    }

    /**
     * {@inheritDoc}
     * ALWAYS RETURNS SELF
     */
    public Location newInstance(Location loc) {
        return loc;
    }

    /**
     * {@inheritDoc}
     * Recursively translates all members of this location.
     */
    public Location translate(int dist) {
        if (this.members.isEmpty()) {
            return this;
        }
        List newmembers = new ArrayList();
        for (Iterator i = this.members.iterator(); i.hasNext();) {
            RichLocation rl = (RichLocation) i.next();
            newmembers.add(rl.translate(dist));
        }
        return new CompoundRichLocation(this.getTerm(), newmembers);
    }

    /**
     * {@inheritDoc}
     * Recursively applies this call to all members. If passed a Location 
     * which is not a RichLocation, it converts it first using 
     * RichLocation.Tools.enrich().
     * @see RichLocation.Tools
     * @return true if an only if one of the members of this <code>Location</code>
     * wholey contains <code>l</code>.
     */
    public boolean contains(Location l) {
        if (!(l instanceof RichLocation)) {
            l = RichLocation.Tools.enrich(l);
        }
        if (l instanceof EmptyRichLocation) {
            return l.contains(this); // let them do the hard work!
        } else {
            RichLocation rl = (RichLocation) l;
            if (rl instanceof CompoundRichLocation) {
                CompoundRichLocation crl = (CompoundRichLocation) rl;
                Map matches = new HashMap();
                for (Iterator i = crl.members.iterator(); i.hasNext();) {
                    matches.put(i.next(), Boolean.FALSE);
                }
                for (Iterator i = this.members.iterator(); i.hasNext();) {
                    RichLocation member = (RichLocation) i.next();
                    for (Iterator j = matches.entrySet().iterator(); j.hasNext();) {
                        Map.Entry entry = (Map.Entry)j.next();
                        if (entry.getValue().equals(Boolean.TRUE)) {
                            continue;
                        }
                        RichLocation match = (RichLocation) entry.getKey();
                        if (member.contains(match)) {
                            entry.setValue(Boolean.TRUE);
                        }
                    }
                }
                for (Iterator i = matches.values().iterator(); i.hasNext(); ) {
                    if (i.next().equals(Boolean.FALSE)) {
                        return false;
                    }
                }
                return true;
            } else {
                for (Iterator i = this.members.iterator(); i.hasNext();) {
                    if (((RichLocation) i.next()).contains(rl)) {
                        return true;
                    }
                }
            }
            return false;
        }
    }

    /**
     * {@inheritDoc}
     * Recursively applies this call to all members.
     * @return true if and only if at least on of the members overlaps <code>l</code>
     */
    public boolean overlaps(Location l) {
        for (Iterator i = this.members.iterator(); i.hasNext();) {
            RichLocation rl = (RichLocation) i.next();
            if (rl.overlaps(l)) {
                return true;
            }
        }
        return false;
    }

    /**
     * {@inheritDoc} 
     * If passed a Location which is not a RichLocation, it converts it first 
     * using RichLocation.Tools.enrich().
     * The resulting location may or may not be a compound location. If it is
     * a compound location, its contents will be a set of simple locations.
     * Regions that overlap will be merged into a single location.
     * @see RichLocation.Tools
     * @return a <code>CompoundLocation</code> if the components of the union
     * cannot be merged else a <code>SimpleRichLocation</code>
     */
    public Location union(Location l) {
        if (!(l instanceof RichLocation)) {
            l = RichLocation.Tools.enrich(l);
        }
        if (l instanceof EmptyRichLocation) {
            return this;
        } else {
            // Easy - construct a new location based on the members of both
            // ourselves and the location passed as a parameter
            List members = new ArrayList();
            members.addAll(RichLocation.Tools.flatten(this));
            members.addAll(RichLocation.Tools.flatten((RichLocation) l));
            return RichLocation.Tools.construct(RichLocation.Tools.merge(members));
        }
    }

    /**
     * {@inheritDoc}
     * If passed a Location which is not a RichLocation, it converts it first 
     * using RichLocation.Tools.enrich().
     * The resulting location may or may not be a compound location. If it is
     * a compound location, its contents will be a set of simple locations.
     * @return a <code>CompoundLocation</code> if there is more than one region
     * of intersection that cannot be merged. Else a <code>SimpleRichLocation</code>.
     */
    public Location intersection(Location l) {
        if (!(l instanceof RichLocation)) {
            l = RichLocation.Tools.enrich(l);
        }
        if (l instanceof EmptyRichLocation) {
            return l;
        } else if (l instanceof CompoundRichLocation) {
            Collection theirMembers = RichLocation.Tools.flatten((RichLocation) l);
            // For every member of the location passed as a parameter, intersect 
            // with ourselves. Then construct a new location from the
            // results of all the intersections.
            Set results = new TreeSet();
            for (Iterator i = theirMembers.iterator(); i.hasNext();) {
                RichLocation member = (RichLocation) i.next();
                results.add(this.intersection(member));
            }
            return RichLocation.Tools.construct(RichLocation.Tools.merge(results));
        } else {
            // Simple vs. ourselves
            // For every member of ourselves, intersect with the location
            // passed as a parameter. Then construct a new location from the
            // results of all the intersections.
            Set results = new TreeSet();
            for (Iterator i = this.members.iterator(); i.hasNext();) {
                RichLocation member = (RichLocation) i.next();
                results.add(member.intersection(l));
            }
            return RichLocation.Tools.construct(RichLocation.Tools.merge(results));
        }
    }

    /**
     * {@inheritDoc}
     * Recursively applies this call to all members.
     */
    public void setCrossRefResolver(CrossReferenceResolver r) {
        for (Iterator i = this.members.iterator(); i.hasNext();) {
            ((RichLocation) i.next()).setCrossRefResolver(r);
        }
    }

    /**
     * {@inheritDoc}
     * This function concatenates the symbols of all its child locations.
     * <p>
     * The most obvious application of this method to a <code>CompoundRichLocation</code>
     * is the contatenation of the components of a gene with multiple exons.
     */
    public SymbolList symbols(SymbolList seq) {
        if (seq == null) {
            throw new IllegalArgumentException("Sequence cannot be null");
        }

        List res = new ArrayList();
        for (Iterator i = this.members.iterator(); i.hasNext();) {
            RichLocation l = (RichLocation) i.next();
            res.addAll(l.symbols(seq).toList());
        }

        try {
            return new SimpleSymbolList(seq.getAlphabet(), res);
        } catch (IllegalSymbolException ex) {
            throw new RuntimeException("Could not build compound sequence string", ex);
        }
    }

    /**
     * {@inheritDoc}
     */
    public void setTerm(ComparableTerm term) throws ChangeVetoException {
        if (term == null) {
            throw new ChangeVetoException("Cannot set term to null");
        }
        super.setTerm(term);
    }

    /**
     * {@inheritDoc}
     */
    public int hashCode() {
        int code = 17;
        code = 31 * code + this.getTerm().hashCode();
        for (Iterator i = this.members.iterator(); i.hasNext();) {
            code = 31 * i.next().hashCode();
        }
        return code;
    }

    /**
     * {@inheritDoc}
     * Compound locations are only equal to other Locations if all their
     * components are equal.
     */
    public boolean equals(Object o) {
        if (o == this) {
            return true;
        }
        if (!(o instanceof Location)) {
            return false;
        }
        Location them = (Location) o;

        if (them.isContiguous()) {
            return false;
        } //because this is not!

        // ok - both compound. The blocks returned from blockIterator should each be
        // equivalent.
        Iterator i1 = this.blockIterator();
        Iterator i2 = them.blockIterator();

        // while there are more pairs to check...
        while (i1.hasNext() && i2.hasNext()) {
            // check that this pair is equivalent
            Location l1 = (Location) i1.next();
            Location l2 = (Location) i2.next();

            if (!(l1.equals(l2))) // not equivalent blocks so not equal
            {
                return false;
            }
        }
        if (i1.hasNext() || i2.hasNext()) {
            // One of the locations had more blocks than the other
            return false;
        }
        // Same number of blocks, all equivalent. Must be equal.
        return true;
    }

    /**
     * {@inheritDoc}
     * 
     */
    public int compareTo(Object o) {
        Location fo = (Location) o;
        if (this.equals(fo)) {
            return 0;
        } else {
            return this.getMin() - fo.getMin();
        }
    }

    /**
     * {@inheritDoc}
     * Form: "term:[location,location...]"
     */
    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(this.getTerm());
        sb.append(":[");
        for (Iterator i = this.blockIterator(); i.hasNext();) {
            sb.append(i.next());
            if (i.hasNext()) {
                sb.append(",");
            }
        }
        sb.append("]");
        return sb.toString();
    }
}

