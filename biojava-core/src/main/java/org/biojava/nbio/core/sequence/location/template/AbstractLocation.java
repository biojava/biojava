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
 * Created on 01-21-2010
 */
package org.biojava.nbio.core.sequence.location.template;

import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.storage.JoiningSequenceReader;
import org.biojava.nbio.core.sequence.template.ComplementCompound;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.sequence.views.ComplementSequenceView;
import org.biojava.nbio.core.sequence.views.ReversedSequenceView;
import org.biojava.nbio.core.util.Hashcoder;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import static java.lang.String.format;
import static org.biojava.nbio.core.util.Equals.classEqual;
import static org.biojava.nbio.core.util.Equals.equal;

/**
 * Base abstraction of a location which encodes for the majority of important
 * features about a location such as the start, end and strand
 *
 * @author ayates
 * @author Paolo Pavan
 */
public abstract class AbstractLocation implements Serializable, Location {

	private static final long serialVersionUID = 1L;

	//TODO Need to have the Sequence lookup resolver here; see the next one as well
	//TODO Need a way of late binding of start/stop

	private Point start;
	private Point end;
	private Strand strand;
	private List<Location> subLocations;
	private boolean circular;
	private boolean betweenCompounds;
	private AccessionID accession;

	private boolean partialOn5prime = false;
	private boolean partialOn3prime = false;



	protected AbstractLocation() {
		super();
	}

	/**
	 * Default constructor
	 *
	 * @param start start of the location
	 * @param end end of the location
	 * @param strand strand it is located on
	 * @param circular Boolean which says if the current location was circular
	 * or not
	 * @param betweenCompounds Indicates the location lies at the position between
	 * a pair of bases; means the bases must be next to each other (and
	 * therefore cannot be complex)
	 * @param subLocations Sub locations which composes this location
	 */
	public AbstractLocation(Point start, Point end, Strand strand,
			boolean circular, boolean betweenCompounds,
			List<Location> subLocations) {
		this(start, end, strand, circular, betweenCompounds, null, subLocations);
	}

	/**
	 * Default constructor
	 *
	 * @param start start of the location
	 * @param end end of the location
	 * @param strand strand it is located on
	 * @param circular Boolean which says if the current location was circular
	 * or not
	 * @param betweenCompounds Indicates the location lies at the position between
	 * a pair of bases; means the bases must be next to each other (and
	 * therefore cannot be complex)
	 * @param accession The accession ID to link this location to
	 * @param subLocations Sub locations which composes this location
	 */
	public AbstractLocation(Point start, Point end, Strand strand,
			boolean circular, boolean betweenCompounds, AccessionID accession,
			List<Location> subLocations) {
		this.start = start;
		this.end = end;
		this.strand = strand;
		this.circular = circular;
		this.betweenCompounds = betweenCompounds;
		this.accession = accession;
		this.subLocations = subLocations==null? null : Collections.unmodifiableList(subLocations);
		assertLocation();
	}

	protected void assertLocation() {
		if (isCircular() && !isComplex()) {
			throw new IllegalStateException("Cannot have a circular "
					+ "location which is not complex");
		}

		int st = getStart().getPosition();
		int e = getEnd().getPosition();

		if (st > e) {
			throw new IllegalStateException(
					String.format("Start (%d) is greater than end (%d); "
					+ "this is an incorrect format",
					st, e));
		}

		if(isBetweenCompounds() && isComplex()) {
			throw new IllegalStateException("Cannot have a complex location "
					+ "which is located between a pair of compounds");
		}

		if(isBetweenCompounds() && (st + 1) != e) {
			throw new IllegalStateException(
					String.format("Start (%d) is not next to end (%d)", st, e));
		}

	}


	@Override
	public Point getEnd() {
		return end;
	}


	@Override
	public Point getStart() {
		return start;
	}


	@Override
	public int getLength() {
		return (getEnd().getPosition() - getStart().getPosition()) + 1;
	}


	@Override
	public Strand getStrand() {
		return strand;
	}


	@Override
	public List<Location> getSubLocations() {
		if(subLocations == null) {
			return Collections.emptyList();
		}
		return subLocations;
	}


	@Override
	public boolean isComplex() {
		return !getSubLocations().isEmpty();
	}


	@Override
	public AccessionID getAccession() {
		return accession;
	}

	public boolean isPartialOn5prime() {
		return partialOn5prime;
	}

	public void setPartialOn5prime(boolean partialOn5prime) {
		this.partialOn5prime = partialOn5prime;
	}

	public boolean isPartialOn3prime() {
		return partialOn3prime;
	}

	public void setPartialOn3prime(boolean partialOn3prime) {
		this.partialOn3prime = partialOn3prime;
	}

	public boolean isPartial() {
		return partialOn5prime || partialOn3prime;
	}

	/**
	 * Iterates through all known sub-locations for this location but does
	 * not descend
	 */

	@Override
	public Iterator<Location> iterator() {
		List<Location> list;
		if(isComplex()) {
			list = getSubLocations();
		}
		else {
			list = new ArrayList<>();
			list.add(this);
		}
		return list.iterator();
	}

	/**
	 * Returns the normalised list of sub locations i.e. only those locations
	 * which do not have a sub location. Useful for when you need to get
	 * the exact elements of a location back for sub sequences.
	 */

	@Override
	public List<Location> getRelevantSubLocations() {
		return getAllSubLocations(this);
	}

	/**
	 * Here to allow for recursion
	 */
	private List<Location> getAllSubLocations(Location location) {
		List<Location> flatSubLocations = new ArrayList<>();
		for (Location l : location.getSubLocations()) {
			if (l.isComplex()) {
				flatSubLocations.addAll(getAllSubLocations(l));
			}
			else {
				flatSubLocations.add(l);
			}
		}
		return flatSubLocations;
	}


	@Override
	public boolean equals(Object obj) {
		boolean equals = false;
		if (classEqual(this, obj)) {
			AbstractLocation l = (AbstractLocation) obj;
			equals = (equal(getStart(), l.getStart())
					&& equal(getEnd(), l.getEnd())
					&& equal(getStrand(), l.getStrand())
					&& equal(isCircular(), l.isCircular())
					&& equal(isBetweenCompounds(), l.isBetweenCompounds())
					&& equal(getSubLocations(), l.getSubLocations())
					&& equal(getAccession(), l.getAccession()));
		}
		return equals;
	}


	@Override
	public int hashCode() {
		int r = Hashcoder.SEED;
		r = Hashcoder.hash(r, getStart());
		r = Hashcoder.hash(r, getEnd());
		r = Hashcoder.hash(r, getStrand());
		r = Hashcoder.hash(r, isCircular());
		r = Hashcoder.hash(r, isBetweenCompounds());
		r = Hashcoder.hash(r, getSubLocations());
		r = Hashcoder.hash(r, getAccession());
		return r;
	}


	@Override
	public boolean isCircular() {
		return circular;
	}


	@Override
	public boolean isBetweenCompounds() {
		return betweenCompounds;
	}

	//TODO Support the accession based lookup system; maybe still require a different impl?

	/**
	 * If circular this will return the sequence represented by the sub
	 * locations joined. If not circular then we get the Sequence for the
	 * outer-bounds defined by this location.
	 */

	@Override
	public <C extends Compound> Sequence<C> getSubSequence(Sequence<C> sequence) {
		if(isCircular()) {
			List<Sequence<C>> sequences = new ArrayList<>();
			for(Location l: this) {
				sequences.add(l.getSubSequence(sequence));
			}
			return new JoiningSequenceReader<>(sequence.getCompoundSet(), sequences);
		}
		return reverseSequence(sequence.getSubSequence(
				getStart().getPosition(), getEnd().getPosition()));
	}

	/**
	 *
	 */

	@Override
	public <C extends Compound> Sequence<C> getRelevantSubSequence(Sequence<C> sequence) {
		List<Sequence<C>> sequences = new ArrayList<>();
		for(Location l: getRelevantSubLocations()) {
			sequences.add(l.getSubSequence(sequence));
		}
		return new JoiningSequenceReader<>(sequence.getCompoundSet(), sequences);
	}

	/**
	 * Reverses and (if possible) complements the Sequence so as to represent
	 * the reverse strand (if one exists). Also does checking to see if the
	 * location we are on is on the reverse strand or not.
	 */
	@SuppressWarnings("unchecked")
	protected <C extends Compound> Sequence<C> reverseSequence(Sequence<C> sequence) {
		if(getStrand() != Strand.NEGATIVE) {
			return sequence;
		}

		Sequence<C> reversed = new ReversedSequenceView<>(sequence);
		// "safe" operation as we have tried to check this
		if(canComplement(sequence)) {
			Sequence<ComplementCompound> casted = (Sequence<ComplementCompound>) reversed;
			ComplementSequenceView<ComplementCompound> complement =
					new ComplementSequenceView<>(casted);
			return (Sequence<C>)complement;
		}
		return reversed;
	}

	/**
	 * Uses the Sequence's CompoundSet to decide if a compound can
	 * be assgined to ComplementCompound meaning it can complement
	 */
	protected <C extends Compound> boolean canComplement(Sequence<C> sequence) {
		CompoundSet<C> compoundSet = sequence.getCompoundSet();
		Compound c = compoundSet.getAllCompounds().iterator().next();
		return ComplementCompound.class.isAssignableFrom(c.getClass());
	}


	@Override
	public String toString() {
		String circ = (isCircular()) ? " - circular" : "";
		String between = (isBetweenCompounds()) ? "^" : "..";
		return format("%d%s%d(%s%s)", getStart().getPosition(), between, getEnd().getPosition(),
				getStrand().getStringRepresentation(), circ);
	}

	protected void setCircular(boolean circular) {
		this.circular = circular;
	}

	protected void setEnd(Point end) {
		this.end = end;
	}

	protected void setStart(Point start) {
		this.start = start;
	}

	public void setStrand(Strand strand) {
		this.strand = strand;
	}

	public void setBetweenCompounds(boolean betweenCompounds) {
		this.betweenCompounds = betweenCompounds;
	}

	public void setSubLocations(List<Location> subLocations) {
		this.subLocations = subLocations;
	}

	public void setAccession(AccessionID accession) {
		this.accession = accession;
	}
}
