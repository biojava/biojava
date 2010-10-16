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
import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.RNATools;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolListViews;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.CrossRef;
import org.biojavax.CrossReferenceResolver;
import org.biojavax.Namespace;
import org.biojavax.RichAnnotation;
import org.biojavax.RichObjectFactory;
import org.biojavax.SimpleRichAnnotation;
import org.biojavax.bio.seq.io.GenbankLocationParser;
import org.biojavax.ontology.ComparableTerm;

/**
 * A simple implementation of RichLocation.
 * @author Richard Holland
 * @author Mark Schreiber
 * @author George Waldon
 * @since 1.5
 */
public class SimpleRichLocation extends AbstractChangeable implements RichLocation {
    
    private CrossRef crossRef;
    private RichAnnotation notes = new SimpleRichAnnotation();
    protected ComparableTerm term;
    private Position min;
    private Position max;
    private PositionResolver pr = RichObjectFactory.getDefaultPositionResolver();
    private CrossReferenceResolver crr = RichObjectFactory.getDefaultCrossReferenceResolver();
    private Strand strand;
    private int rank;
    protected int circularLength = 0;
    private RichFeature feature;
    
    /**
     * Creates a new instance of SimpleRichSequenceLocation that points to a
     * single position on the positive strand.
     * @param pos the location position (a point).
     * @param rank Rank of location.
     */
    public SimpleRichLocation(Position pos, int rank) {
        this(pos,pos,rank,Strand.POSITIVE_STRAND);
    }
    
    /**
     * Creates a new instance of SimpleRichSequenceLocation that points to a
     * single position.
     * @param pos the location position (a point).
     * @param rank Rank of location.
     * @param strand The strand of the location
     */
    public SimpleRichLocation(Position pos, int rank, Strand strand) {
        this(pos,pos,rank,strand,null);
    }
    
    /**
     * Creates a new instance of SimpleRichSequenceLocation that points to a
     * single position on another sequence.
     * @param pos the location position (a point).
     * @param rank Rank of location.
     * @param strand the strand of the location
     * @param crossRef a cross reference to another object (null for parent sequence)
     */
    public SimpleRichLocation(Position pos, int rank, Strand strand, CrossRef crossRef) {
        this(pos,pos,rank,strand,crossRef);
    }
    
    /**
     * Creates a new instance of SimpleRichSequenceLocation that points to a
     * range position on the positive strand.
     * @param min the minimum bound of the location
     * @param max the maximum bound of the location
     * @param rank Rank of location.
     *
     */
    public SimpleRichLocation(Position min, Position max, int rank) {
        this(min,max,rank,Strand.POSITIVE_STRAND);
    }
    
    /**
     * Creates a new instance of SimpleRichSequenceLocation that points to a
     * range position.
     * @param min the minimum bound of the location
     * @param max the maximum bound of the location
     * @param rank Rank of location.
     * @param strand the strand of the location
     */
    public SimpleRichLocation(Position min, Position max, int rank, Strand strand) {
        this(min,max,rank,strand,null);
    }
    
    /**
     * Creates a new instance of SimpleRichSequenceLocation that points to a
     * range position on another sequence.
     * @param min the minimum bound of the location
     * @param max the maximum bound of the location
     * @param rank Rank of location.
     * @param strand the strand of the location
     * @param crossRef a cross reference to another object (null for parent sequence)
     */
    public SimpleRichLocation(Position min, Position max, int rank, Strand strand, CrossRef crossRef) {
        this.min = min;
        this.max = max;
        this.rank = rank;
        this.strand = strand;
        this.crossRef = crossRef;
        this.feature = null;
    }
    
    // Hibernate requirement - not for public use.
    protected SimpleRichLocation() {}
    
    /**
     * {@inheritDoc}
     */
    public void sort() {}
    
    /**
     * {@inheritDoc}
     */
    public RichFeature getFeature() { return this.feature; }
    
    /**
     * {@inheritDoc}
     */
    public void setFeature(RichFeature feature) throws ChangeVetoException {
        if(!this.hasListeners(RichLocation.FEATURE)) {
            this.feature = feature;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    RichLocation.FEATURE,
                    feature,
                    this.feature
                    );
            ChangeSupport cs = this.getChangeSupport(RichLocation.FEATURE);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.feature = feature;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public CrossRef getCrossRef() { return this.crossRef; }
    
    
    // Hibernate requirement - not for public use.
    protected void setCrossRef(CrossRef crossRef) { this.crossRef = crossRef; }
    
    /**
     * {@inheritDoc}
     */
    public Annotation getAnnotation() { return getRichAnnotation(); }
    
    /**
     * {@inheritDoc}
     */
    public RichAnnotation getRichAnnotation() { return this.notes; }
    
    /**
     * {@inheritDoc}
     * <b>Warning</b> this method gives access to the original
     * Collection not a copy. This is required by Hibernate. If you
     * modify the object directly the behaviour may be unpredictable.
     */
    public Set getNoteSet() { return this.notes.getNoteSet(); }
    
    /**
     * {@inheritDoc}
     * <b>Warning</b> this method gives access to the original
     * Collection not a copy. This is required by Hibernate. If you
     * modify the object directly the behaviour may be unpredictable.
     */
    public void setNoteSet(Set notes) throws ChangeVetoException { this.notes.setNoteSet(notes); }
    
    /**
     * {@inheritDoc}
     */
    public ComparableTerm getTerm() { return this.term; }
    
    /**
     * {@inheritDoc}
     */
    public void setTerm(ComparableTerm term) throws ChangeVetoException {
        if(!this.hasListeners(RichLocation.TERM)) {
            this.term = term;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    RichLocation.TERM,
                    term,
                    this.term
                    );
            ChangeSupport cs = this.getChangeSupport(RichLocation.TERM);
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
    public int getCircularLength() { return this.circularLength; }
    
    /**
     * {@inheritDoc}
     */
    public void setCircularLength(int circularLength) throws ChangeVetoException {
        if(!this.hasListeners(RichLocation.CIRCULAR)) {
            this.circularLength = circularLength;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    RichLocation.CIRCULAR,
                    new Integer(circularLength),
                    new Integer(this.circularLength)
                    );
            ChangeSupport cs = this.getChangeSupport(RichLocation.CIRCULAR);
            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                this.circularLength = circularLength;
                cs.firePostChangeEvent(ce);
            }
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public Strand getStrand() { return this.strand; }
    
    // Hibernate requirement - not for public use.
    protected void setStrand(Strand strand) { this.strand = strand; }
    
    // Hibernate requirement - not for public use.
    int getStrandNum() { return this.strand.intValue(); }
    
    // Hibernate requirement - not for public use.
    void setStrandNum(int token) { this.strand = Strand.forValue(token); }
    
    /**
     * {@inheritDoc}
     */
    public int getRank() { return this.rank; }
    
    /**
     * {@inheritDoc}
     */
    public void setRank(int rank) throws ChangeVetoException {
        if(!this.hasListeners(RichLocation.RANK)) {
            this.rank = rank;
        } else {
            ChangeEvent ce = new ChangeEvent(
                    this,
                    RichLocation.RANK,
                    new Integer(rank),
                    new Integer(this.rank)
                    );
            ChangeSupport cs = this.getChangeSupport(RichLocation.RANK);
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
    public int getMax() {
        if (this.max.equals(this.min)) return this.getMin(); // because resolver might resolve differently
        else return this.pr.getMax(this.max);
    }
    
    // Hibernate requirement - not for public use.
    void setMax(int max) {  this.max = new SimplePosition(false,false,max); }
    
    /**
     * {@inheritDoc}
     */
    public int getMin() { return this.pr.getMin(this.min); }
    
    // Hibernate requirement - not for public use.
    void setMin(int min) {  this.min = new SimplePosition(false,false,min); }
    
    /**
     * {@inheritDoc}
     */
    public Position getMinPosition() { return this.min; }
    
    // Hibernate requirement - not for public use.
    protected void setMinPosition(Position min) {  this.min = min; }
    
    /**
     * {@inheritDoc}
     */
    public Position getMaxPosition() { return this.max; }
    
    // Hibernate requirement - not for public use.
    protected void setMaxPosition(Position max) {  this.max = max; }
    
    /**
     * {@inheritDoc}
     */
    public void setPositionResolver(PositionResolver p) { this.pr = p; }
    
    /**
     * {@inheritDoc}
     */
    public Iterator blockIterator() { return Collections.singleton(this).iterator(); }
    
    /**
     * {@inheritDoc}
     */
    public boolean isContiguous() { return true; }
    
    /**
     * {@inheritDoc}
     */
    public boolean contains(int p) {
        int modStart = this.getMin();
        int modEnd = this.getMax();
        if (this.circularLength>0) {
            // Modulate the point to fall inside our sequence
            p = RichLocation.Tools.modulateCircularIndex(p, this.circularLength);
            // Modulate our location to the earliest possible point in our sequence
            int[] ourModParts = RichLocation.Tools.modulateCircularLocation(modStart, modEnd, this.circularLength);
            modStart = ourModParts[0];
            modEnd = ourModParts[1];
            // If we wrap and the point occurs before us, increment point to fall in correct range
            if (modEnd>this.circularLength && p<modStart) p+=this.circularLength;
        }
        return (p>=modStart && p<=modEnd);
    }
    
    /**
     * {@inheritDoc}
     */
    public Location getDecorator(Class decoratorClass) { return null; }
    
    /**
     * {@inheritDoc}
     */
    public Location newInstance(Location loc) { return loc; }
    
    /**
     * {@inheritDoc}
     */
    public Location translate(int dist) {
        return new SimpleRichLocation(this.min.translate(dist),this.max.translate(dist),0,this.strand,this.crossRef);
    }
    
    /**
     * {@inheritDoc}
     * A location contains another location if it overlaps it, and the coordinates
     * enclose those of the other location at both ends, and they fall on
     * the same strand.
     */
    public boolean contains(Location l) {
        if (!(l instanceof RichLocation)) l = RichLocation.Tools.enrich(l);
        if (l instanceof EmptyRichLocation) return false;
        else {
            RichLocation rl = (RichLocation)l;
            // Simple vs. simple
            if (!this.overlaps(rl)) return false; // No overlap = no possible contain
            if (!this.getStrand().equals(rl.getStrand())) return false; // Diff strand = not contained
            if (this.circularLength>0) {
                // Modulate to shortest possible equivalent region
                int parts[] = RichLocation.Tools.modulateCircularLocationPair(this,rl,this.circularLength);
                int ourModStart = parts[0];
                int ourModEnd = parts[1];
                int theirModStart = parts[2];
                int theirModEnd = parts[3];
                return (ourModStart <= theirModStart && ourModEnd >= theirModEnd);
            } else {
                return (this.getMin() <= rl.getMin() && this.getMax() >= rl.getMax());
            }
        }
    }
    
    /**
     * {@inheritDoc}
     * A location overlaps another location if it is on the same sequence, and
     * the coordinates overlap, and both are of the same circularity.
     */
    public boolean overlaps(Location l) {
        if (!(l instanceof RichLocation)) l = RichLocation.Tools.enrich(l);
        if (l instanceof EmptyRichLocation) return false;
        else if (l instanceof CompoundRichLocation) return l.overlaps(this);
        else {
            // Simple vs. simple
            RichLocation rl = (RichLocation)l;
            // Mismatch of crossref is no overlap.
            if (rl.getCrossRef()!=null || this.crossRef!=null) {
                if (rl.getCrossRef()!=null && this.crossRef!=null) {
                    if (!this.crossRef.equals(rl.getCrossRef())) return false;
                } else return false;
            }
            if (this.circularLength!=rl.getCircularLength()) return false; // Diff circularLength location sizes = not overlapping
            // Modulate our start/end to shortest possible equivalent region
            if (this.circularLength>0) {
                int parts[] = RichLocation.Tools.modulateCircularLocationPair(this,rl,this.circularLength);
                int ourModStart = parts[0];
                int ourModEnd = parts[1];
                int theirModStart = parts[2];
                int theirModEnd = parts[3];
                return (ourModStart<=theirModEnd && ourModEnd>=theirModStart);
            } else {
                return (this.getMin()<=rl.getMax() && this.getMax()>=rl.getMin());
            }
        }
    }
    
    /**
     * {@inheritDoc}
     * A merged SimpleRichLocation is returned if the two locations overlap and are on
     * the same strand. Otherwise, a CompoundRichLocation is returned containing
     * the two locations as members.
     */
    public Location union(Location l) {
        if (!(l instanceof RichLocation)) l = RichLocation.Tools.enrich(l);
        if (l instanceof EmptyRichLocation) return this;
        else if (l instanceof CompoundRichLocation) return l.union(this);
        else {
            // Simple vs. simple
            RichLocation rl = (RichLocation)l;
            if (this.overlaps(rl) && this.getStrand().equals(rl.getStrand())) {
                // We can do the one-v-one overlapping same-strand union
                if (this.circularLength>0) {
                    // Union of Overlapping circular locations
                    // Modulate our start/end to shortest possible equivalent region
                    int parts[] = RichLocation.Tools.modulateCircularLocationPair(this,rl,this.circularLength);
                    int ourModStart = parts[0];
                    int ourModEnd = parts[1];
                    int theirModStart = parts[2];
                    int theirModEnd = parts[3];
                    // Now we can select the minimum and maximum positions using the modded locations
                    Position startPos = (ourModStart<theirModStart)?this.min:rl.getMinPosition();
                    Position endPos = (ourModEnd>theirModEnd)?this.max:rl.getMaxPosition();
                    return new SimpleRichLocation(startPos,endPos,0,this.strand,this.crossRef);
                } else {
                    // Union of Overlapping non-circular locations
                    return new SimpleRichLocation(this.posmin(this.min,rl.getMinPosition()),this.posmax(this.max,rl.getMaxPosition()),0,this.strand,this.crossRef);
                }
            }
            // We can do the one-v-one non-overlapping or different-strand union too
            Collection members = new ArrayList();
            members.add(this);
            members.add(l);
            if (RichLocation.Tools.isMultiSource(members)) return new MultiSourceCompoundRichLocation(members);
            else return new CompoundRichLocation(members);
        }
    }
    
    /**
     * {@inheritDoc}
     * If the locations overlap and are on the same strand, the intersection
     * is returned. If they overlap but are on different strands, a CompoundRichLocation
     * of the overlapping portions is returned. If they do not overlap, the empty
     * sequence is returned.
     */
    public Location intersection(Location l) {
        if (!(l instanceof RichLocation)) l = RichLocation.Tools.enrich(l);
        if (l instanceof EmptyRichLocation) return l;
        else if (l instanceof CompoundRichLocation) return l.intersection(this);
        else {
            RichLocation rl = (RichLocation)l;
            if (this.overlaps(l)) {
                if (this.getStrand().equals(rl.getStrand())) {
                    // We can do the one-v-one same-strand overlapping intersection here
                    if (this.circularLength>0) {
                        // Modulate our start/end to shortest possible equivalent region
                        int parts[] = RichLocation.Tools.modulateCircularLocationPair(this,rl,this.circularLength);
                        int ourModStart = parts[0];
                        int ourModEnd = parts[1];
                        int theirModStart = parts[2];
                        int theirModEnd = parts[3];
                        // Now we can select the minimum and maximum positions using the modded locations
                        Position startPos = (ourModStart>theirModStart)?this.min:rl.getMinPosition();
                        Position endPos = (ourModEnd<theirModEnd)?this.max:rl.getMaxPosition();
                        return new SimpleRichLocation(startPos,endPos,0,this.strand,this.crossRef);
                    } else {
                        return new SimpleRichLocation(this.posmax(this.min,rl.getMinPosition()),this.posmin(this.max,rl.getMaxPosition()),0,this.strand,this.crossRef);
                    }
                }
                // We can do the one-v-one different-strand overlapping intersection here
                Collection members = new ArrayList();
                members.add(new SimpleRichLocation(this.posmax(this.min,rl.getMinPosition()),this.posmin(this.max,rl.getMaxPosition()),0,this.strand,this.crossRef));
                members.add(new SimpleRichLocation(this.posmax(this.min,rl.getMinPosition()),this.posmin(this.max,rl.getMaxPosition()),0,rl.getStrand(),this.crossRef));
                if (RichLocation.Tools.isMultiSource(members)) return new MultiSourceCompoundRichLocation(members);
                else return new CompoundRichLocation(members);
            } else {
                // We can do the one-v-one non-overlapping intersection here
                return RichLocation.EMPTY_LOCATION;
            }
        }
    }
    
    // calculates the smaller of the two positions, based on their resolver output
    protected Position posmin(Position a, Position b) {
        int ar = this.pr.getMin(a);
        int br = this.pr.getMin(b);
        if (ar<=br) return a;
        else return b;
    }
    
    // calculates the smaller of the two positions, based on their resolver output
    protected Position posmax(Position a, Position b) {
        int ar = this.pr.getMax(a);
        int br = this.pr.getMax(b);
        if (ar>br) return a;
        else return b;
    }
    
    /**
     * {@inheritDoc}
     */
    public void setCrossRefResolver(CrossReferenceResolver r) {
        if (r==null) throw new IllegalArgumentException("Resolver cannot be null");
        this.crr = r;
    }
    
    /**
     * {@inheritDoc}
     * If the location is circular but the sequence is not, or they are both
     * circular but of different circular lengths, an exception is thrown.
     * The symbol list passed in is the sequence used to obtain symbols
     * if the cross reference for this location has not been set. If the cross
     * reference has been set, then the symbol list passed in is only used
     * if it has the same accession, namespace and version as the cross
     * reference on this location. Otherwise, the cross referenced symbol list
     * is looked up and used instead.
     */
    public SymbolList symbols(SymbolList seq) {
        if (seq==null) throw new IllegalArgumentException("Sequence cannot be null");
        
        if (seq instanceof RichSequence) {
            RichSequence rs = (RichSequence)seq;
            if (this.getCircularLength()>0) {
                if (!rs.getCircular()) throw new IllegalArgumentException("Attempt to apply circular location to non-circular sequence");
                if (rs.length()==this.getCircularLength()) throw new IllegalArgumentException("Attempt to apply circular location to circular sequence of different length");
            }
        }
        
        // Resolve cross-references to remote sequences
        if (this.getCrossRef()!=null) {
            CrossRef cr = this.getCrossRef();
            if (seq instanceof RichSequence) {
                RichSequence rs = (RichSequence)seq;
                String accession = rs.getAccession();
                Namespace ns = rs.getNamespace();
                String raccession = cr.getAccession();
                String rnamespace = cr.getDbname();
                if (!(accession.equals(raccession) && ns.getName().equals(rnamespace))) {
                    // It really is remote - the xref doesn't point to the sequence we just got passed
                    seq = this.crr.getRemoteSymbolList(cr, seq.getAlphabet());
                }
            } else {
                // It's assumed to be remote because we can't tell what the sequence we were passed really is
                seq = this.crr.getRemoteSymbolList(this.getCrossRef(), seq.getAlphabet());
            }
        }
        
        
        // Carry on as before
        seq = seq.subList(this.getMin(),this.getMax());
        
        try {
            if (this.strand==Strand.NEGATIVE_STRAND) {
                Alphabet a = seq.getAlphabet();
                if (a==AlphabetManager.alphabetForName("DNA")) {
                    seq = DNATools.reverseComplement(seq);
                } else if (a==AlphabetManager.alphabetForName("RNA")) {
                    seq = RNATools.reverseComplement(seq);
                } else {
                    seq = SymbolListViews.reverse(seq);// no complement as no such thing
                }
            }
        } catch (IllegalAlphabetException e) {
            IllegalArgumentException ex =
                    new IllegalArgumentException("Could not understand alphabet of passed sequence");
            ex.initCause(e);
            throw ex;
        }
        
        return seq;
    }
    
    /**
     * {@inheritDoc}
     */
    public int hashCode() {
        int code = 17;
        // Hibernate comparison - we haven't been populated yet
        if (this.strand==null) return code;
        // Normal comparison
        if (this.term!=null) code = 31*code + this.term.hashCode();
        code = 31*code + this.getMin();
        code = 31*code + this.getMax();
        code = 31*code + this.strand.hashCode();
        code = 31*code + this.rank;
        if (this.crossRef!=null) code = 31*code + this.crossRef.hashCode();
        return code;
    }
    
    /**
     * {@inheritDoc}
     * Locations are equal if their term, min, max, strand, and crossref are
     * the same, and if their rank is the same too.
     */
    public boolean equals(Object o) {
        if (! (o instanceof RichLocation)) return false;
        // Hibernate comparison - we haven't been populated yet
        if (this.strand==null) return false;
        // Normal comparison
        RichLocation fo = (RichLocation) o;
        if (this.term!=null || fo.getTerm()!=null) {
            if (this.term!=null && fo.getTerm()!=null) {
                if (!this.term.equals(fo.getTerm())) return false;
            } else return false;
        }
        if (this.getMin()!=fo.getMin()) return false;
        if (this.getMax()!=fo.getMax()) return false;
        if (!this.strand.equals(fo.getStrand())) return false;
        if (this.crossRef!=null || fo.getCrossRef()!=null) {
            if (this.crossRef!=null && fo.getCrossRef()!=null) {
                if(!this.crossRef.equals(fo.getCrossRef())) return false;
            } else return false;
        }
        return this.rank==fo.getRank();
    }
    
    /**
     * {@inheritDoc}
     * Locations are sorted first by rank, then crossref, then
     * strand, then term, then min, then max.
     */
    public int compareTo(Object o) {
        if (o==this) return 0;
        // Hibernate comparison - we haven't been populated yet
        if (this.strand==null) return -1;
        // Check if we can really compare at all
        if (!(o instanceof RichLocation)) return -1;
        // Normal comparison
        RichLocation fo = (RichLocation) o;
        if (this.rank!=fo.getRank()) return this.rank-fo.getRank();
        if (this.crossRef!=null || fo.getCrossRef()!=null) {
            if (this.crossRef!=null && fo.getCrossRef()!=null) {
                return this.crossRef.compareTo(fo.getCrossRef());
            } else return -1;
        }
        if (!this.strand.equals(fo.getStrand())) return this.strand.compareTo(fo.getStrand());
        if (this.term!=null || fo.getTerm()!=null) {
            if (this.term!=null && fo.getTerm()!=null && !this.term.equals(fo.getTerm())) return this.term.compareTo(fo.getTerm());
            else return -1;
        }
        if (this.getMin()!=fo.getMin()) return this.getMin()-fo.getMin();
        return this.getMax()-fo.getMax();
    }
    
    /**
     * {@inheritDoc}
     * Form: "start..end" or just "point" for point locations
     */
    public String toString() {
        if (this.max.equals(this.min)) {
            return this.min.toString();
        } else {
            return this.min+".."+this.max;
        }
    }
    
    // Internal use only.
    void setLocationText(final String theLocation) throws ParseException {
//    	System.out.println("SimpleRichLocation.setLocationText-theLocation: ["+theLocation+"]");
    	if (theLocation == null) {
    		setMinPosition(RichLocation.EMPTY_LOCATION.getMinPosition());
    		setMaxPosition(RichLocation.EMPTY_LOCATION.getMaxPosition());
    	} else {
        	final RichLocation location = GenbankLocationParser.parseLocation(RichObjectFactory.getDefaultNamespace(), null, theLocation);
    		setMinPosition(location.getMinPosition());
    		setMaxPosition(location.getMaxPosition());
    	}
//    	System.out.println("SimpleRichLocation.setLocationText-this: ["+this+"]");
    }
    
    // Internal use only.
    String getLocationText() {
//    	System.out.println("SimpleRichLocation.getLocationText-returns: ["+GenbankLocationParser.writeLocation(new SimpleRichLocation(getMinPosition(), getMaxPosition(), getRank()))+"], this: ["+this+"]");
    	return  GenbankLocationParser.writeLocation(new SimpleRichLocation(getMinPosition(), getMaxPosition(), getRank()));
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

