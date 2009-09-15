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

import org.biojava.bio.symbol.FuzzyLocation;
import org.biojava.bio.symbol.FuzzyPointLocation;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.MergeLocation;
import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.CrossRef;
import org.biojavax.CrossReferenceResolver;
import org.biojavax.RichAnnotatable;
import org.biojavax.ontology.ComparableTerm;

/**
 * Describes locations, and adds the concepts of circularity, fuzziness,
 * annotations, and cross-references to other databases. Also includes strands.
 * 
 * @author Richard Holland
 * @since 1.5
 */
public interface RichLocation extends Location, RichAnnotatable, Comparable {

	public static final ChangeType NOTE = new ChangeType(
			"This location's notes have changed",
			"org.biojavax.bio.seq.RichLocation", "NOTE");
	public static final ChangeType TERM = new ChangeType(
			"This location's term has changed",
			"org.biojavax.bio.seq.RichLocation", "TERM");
	public static final ChangeType RANK = new ChangeType(
			"This location's rank has changed",
			"org.biojavax.bio.seq.RichLocation", "RANK");
	public static final ChangeType CIRCULAR = new ChangeType(
			"This location's circularity has changed",
			"org.biojavax.bio.seq.RichLocation", "CIRCULAR");
	public static final ChangeType FEATURE = new ChangeType(
			"This location's parent feature has changed",
			"org.biojavax.bio.seq.RichLocation", "FEATURE");

	/**
	 * The empty location matches nothing.
	 */
	public static final RichLocation EMPTY_LOCATION = new EmptyRichLocation();

	/**
	 * Sorts the member locations of a compound location. Does nothing for
	 * non-compound locations. Sorting depends on the compareTo() method of the
	 * member locations - usually they will be sorted by their start position.
	 * This might be useful when used with the location returned by a union or
	 * intersection, or when setting the term of a compound location to
	 * ORDER_TERM.
	 */
	public void sort();

	/**
	 * Retrieves the feature this location is associated with. May be null.
	 * 
	 * @return the feature.
	 */
	public RichFeature getFeature();

	/**
	 * Sets the feature this location is associated with. If null, that's fine,
	 * but you won't be able to persist it to the database until you give it a
	 * not-null value.
	 * 
	 * @param feature
	 *            the feature.
	 */
	public void setFeature(RichFeature feature) throws ChangeVetoException;

	/**
	 * Retrieves the crossref associated with this location.
	 * 
	 * @return the crossref.
	 */
	public CrossRef getCrossRef();

	/**
	 * Retrieves the term associated with this location.
	 * 
	 * @return the term.
	 */
	public ComparableTerm getTerm();

	/**
	 * Sets the term for this location.
	 * 
	 * @param term
	 *            the term this location should adopt.
	 * @throws ChangeVetoException
	 *             in case of error.
	 */
	public void setTerm(ComparableTerm term) throws ChangeVetoException;

	/**
	 * Retrieves the strand associated with this location.
	 * 
	 * @return the strand.
	 */
	public Strand getStrand();

	/**
	 * Retrieves the rank associated with this location.
	 * 
	 * @return the rank.
	 */
	public int getRank();

	/**
	 * Sets the rank for this location.
	 * 
	 * @param rank
	 *            the rank this location should adopt.
	 * @throws ChangeVetoException
	 *             in case of error.
	 */
	public void setRank(int rank) throws ChangeVetoException;

	/**
	 * Retrieves the start position of this location.
	 * 
	 * @return the position.
	 */
	public Position getMinPosition();

	/**
	 * Retrieves the end position of this location.
	 * 
	 * @return the position.
	 */
	public Position getMaxPosition();

	/**
	 * Sets the resolver to use when working out actual base coordinates from
	 * fuzzy positions.
	 * 
	 * @param p
	 *            the position resolver to use.
	 */
	public void setPositionResolver(PositionResolver p);

	/**
	 * Retrieves the circular length of this location. If it is 0, the location
	 * is not circular. If it is not zero, then the number refers to the
	 * wrapping length of the location. eg. 100 would signify that a position of
	 * 112 would actually be a position of 112-100 = 12.
	 * 
	 * @return the position.
	 */
	public int getCircularLength();

	/**
	 * Sets the circular length of this location. If it is 0, the location is
	 * not circular. If it is not zero, then the number refers to the wrapping
	 * length of the location. eg. 100 would signify that a position of 112
	 * would actually be a position of 112-100 = 12.
	 * 
	 * @param sourceSeqLength
	 *            the circular length of this location
	 * @throws ChangeVetoException
	 *             if it doesn't want to change.
	 */
	public void setCircularLength(int sourceSeqLength)
			throws ChangeVetoException;

	/**
	 * Sets the cross ref resolver to use when retrieving remote symbols. If
	 * none is given, then the default from
	 * RichObjectFactory.getDefaultCrossRefResolver() is used.
	 * 
	 * @param r
	 *            the resolver to use.
	 */
	public void setCrossRefResolver(CrossReferenceResolver r);

	/**
	 * This class represents a strand on which a location may lie. Three strands
	 * are defined by default - UNKNOWN, NEGATIVE, and POSITIVE.
	 */
	public static class Strand implements Comparable {

		private String name;
		private int value;

		/**
		 * The positive strand is represented by the symbol '+' and has the
		 * number 1.
		 */
		public static final Strand POSITIVE_STRAND = new Strand("+", 1);

		/**
		 * The negative strand is represented by the symbol '-' and has the
		 * number -1.
		 */
		public static final Strand NEGATIVE_STRAND = new Strand("-", -1);

		/**
		 * The unknown strand is represented by the symbol '?' and has the
		 * number 0.
		 */
		public static final Strand UNKNOWN_STRAND = new Strand("?", 0);

		/**
		 * Returns the strand object that matches the number given. Throws an
		 * exception if it could not recognise the number. Number is usually
		 * 1,-1,0.
		 * 
		 * @param value
		 *            the number of the strand.
		 * @return the strand matching that number.
		 */
		public static Strand forValue(int value) {
			switch (value) {
			case 1:
				return POSITIVE_STRAND;
			case 0:
				return UNKNOWN_STRAND;
			case -1:
				return NEGATIVE_STRAND;
			default:
				throw new IllegalArgumentException("Unknown strand type: "
						+ value);
			}
		}

		/**
		 * Returns the strand object that matches the symbol given. Throws an
		 * exception if it could not recognise the symbol. Symbol is usually
		 * +,-,?.
		 * 
		 * @param name
		 *            the symbol of the strand.
		 * @return the strand matching that symbol.
		 */
		public static Strand forName(String name) {
			if (name.equals("+"))
				return POSITIVE_STRAND;
			else if (name.equals("?"))
				return UNKNOWN_STRAND;
			else if (name.equals("."))
				return UNKNOWN_STRAND;
			else if (name.equals("-"))
				return NEGATIVE_STRAND;
			else
				throw new IllegalArgumentException("Unknown strand type: "
						+ name);
		}

		// creates a strand with the given number and value
		private Strand(String name, int value) {
			this.name = name;
			this.value = value;
		}

		/**
		 * Returns the numeric value of this strand.
		 * 
		 * @return the numeric value.
		 */
		public int intValue() {
			return this.value;
		}

		/**
		 * Returns the string symbol of this strand.
		 * 
		 * @return the string symbol.
		 */
		public String getName() {
			return this.name;
		}

		/**
		 * {@inheritDoc} Form: "symbol" (eg. +,-,?)
		 */
		public String toString() {
			return this.name;
		}

		/**
		 * {@inheritDoc}
		 */
		public int hashCode() {
			int code = 17;
			code = 31 * code + this.name.hashCode();
			code = 31 * code + this.value;
			return code;
		}

		/**
		 * {@inheritDoc} Strands are equal if their numbers and symbols match.
		 */
		public boolean equals(Object o) {
			if (!(o instanceof Strand))
				return false;
			if (o == this)
				return true;
			Strand them = (Strand) o;
			if (!them.toString().equals(this.name))
				return false;
			if (them.intValue() != this.value)
				return false;
			return true;
		}

		/**
		 * {@inheritDoc} Strands are compared first by symbol, then by number.
		 */
		public int compareTo(Object o) {
			Strand fo = (Strand) o;
			if (!this.name.equals(fo.toString()))
				return this.name.compareTo(fo.toString());
			return this.value - fo.intValue();
		}
	}

	/**
	 * Some useful tools for working with Locations.
	 */
	public static class Tools {

		// because we are static, we don't want to get instantiated
		private Tools() {
		}

		/**
		 * Constructs a RichLocation object based on the given collection of
		 * members. It the collection contains a single location, that is
		 * returned. If it contains multiple locations it returns a
		 * CompoundRichLocation covering them all, with the default term
		 * associated. It returns the empty location if the collection was
		 * empty.
		 * 
		 * @param members
		 *            the members to construct a location from.
		 * @return the corresponding RichLocation
		 */
		public static RichLocation construct(Collection<Location> members) {
			if (members.size() == 0)
				return RichLocation.EMPTY_LOCATION;
			else if (members.size() == 1)
				return members.toArray(new SimpleRichLocation[0])[0];
			else if (isMultiSource(members))
				return new MultiSourceCompoundRichLocation(members);
			else
				return new CompoundRichLocation(members);
		}

		/**
		 * Returns false if all the locations in the set are from the same
		 * strand of the same sequence.
		 * 
		 * @param members
		 *            the set of locations to check.
		 * @return true if they are from multiple sources.
		 */
		public static boolean isMultiSource(Collection<Location> members) {
			RichLocation previous = null;
			for (Iterator<Location> i = members.iterator(); i.hasNext();) {
				RichLocation rl = enrich(i.next());
				if (previous == null)
					previous = rl;
				else {
					if (previous.getCircularLength() != rl.getCircularLength())
						return true;
					if ((previous.getCrossRef() == null && rl.getCrossRef() != null)
							|| (previous.getCrossRef() != null && rl
									.getCrossRef() == null)
							|| (previous.getCrossRef() != rl.getCrossRef() && !previous
									.getCrossRef().equals(rl.getCrossRef())))
						return true;
					if ((previous.getStrand() == null && rl.getStrand() != null)
							|| (previous.getStrand() != null && rl.getStrand() == null)
							|| (previous.getStrand() != rl.getStrand() && !previous
									.getStrand().equals(rl.getStrand())))
						return true;
				}
			}
			return false;
		}

		/**
		 * Takes a set of locations and tries to merge all pairs where the union
		 * operation results in a simple rich location, not a compound one.
		 * 
		 * @param members
		 *            the members to merge
		 * @return the resulting merged set, which may have only one location in
		 *         it.
		 */
		public static Collection<Location> merge(Collection<Location> members) {
			// flatten them out first so we don't end up recursing
			List<Location> membersList = new ArrayList<Location>(
					flatten(members));
			// all members are now singles so we can use single vs single union
			// operations
			if (membersList.size() > 1) {
				for (int p = 0; p < (membersList.size() - 1); p++) {
					RichLocation parent = (RichLocation) membersList.get(p);
					for (int c = p + 1; c < membersList.size(); c++) {
						RichLocation child = (RichLocation) membersList.get(c);
						RichLocation union = (RichLocation) parent.union(child);
						// if parent can merge with child
						if (union.isContiguous()) {
							// replace parent with union
							membersList.set(p, union);
							// remove child
							membersList.remove(c);
							// check all children again
							c = p + 1;
						}
					}
				}
			}
			return membersList;
		}

		/**
		 * Takes a location and returns the set of all members. If any members
		 * are compound, it flattens them too.
		 * 
		 * @param location
		 *            the location to flatten
		 * @return the flattened collection of members.
		 */
		public static Collection<Location> flatten(RichLocation location) {
			List<Location> members = new ArrayList<Location>();
			for (Iterator<Location> i = location.blockIterator(); i.hasNext();)
				members.add(i.next());
			return flatten(members);
		}

		/**
		 * Takes a set of locations and returns the set of all members. If any
		 * members are compound, it flattens them too.
		 * 
		 * @param members
		 *            the locations to flatten
		 * @return the flattened collection of members.
		 */
		public static Collection<Location> flatten(Collection<Location> members) {
			List<Location> flattened = new ArrayList<Location>(members);
			for (int i = 0; i < flattened.size(); i++) {
				RichLocation member = (RichLocation) flattened.get(i);
				if (!member.isContiguous()) {
					flattened.remove(i);
					int insertPos = i;
					for (Iterator<Location> j = member.blockIterator(); j
							.hasNext();)
						flattened.add(insertPos++, j.next());
					i--;
				}
			}
			return flattened;
		}

		/**
		 * Takes a start and end position on a circular location of given
		 * length, and shifts them left along the sequence until they sit at the
		 * earliest possible point where they still would represent the same
		 * sequence.
		 * 
		 * @param start
		 *            the start of the circular location
		 * @param end
		 *            the end of the circular location
		 * @param seqLength
		 *            the circular length of the sequence underlying the
		 *            location
		 * @return an integer array where [0] is the translated start and [1]
		 *         the end.
		 */
		public static int[] modulateCircularLocation(int start, int end,
				int seqLength) {
			// Dummy case for non-circular sequences.
			if (seqLength == 0)
				return new int[] { start, end };
			// Move the end to after the start.
			while (end < start)
				end += seqLength;
			// Calculate the length.
			int locationLength = end - start;
			// Move the start back till it can go no further
			while (start >= seqLength)
				start -= seqLength;
			// Move the end back.
			end = start + locationLength;
			// Return results.
			return new int[] { start, end };
		}

		/**
		 * Takes two circular locations of given length, and shifts them left
		 * along the sequence until they sit at the earliest possible point
		 * where they still would represent the same sequence. The end result
		 * ensures that simple overlap calculations will always work on the
		 * coordinates returned.
		 * 
		 * @param a
		 *            the first location to shift
		 * @param b
		 *            the second location to shift
		 * @param seqLength
		 *            the circular length of the sequence underlying the
		 *            location
		 * @return an integer array where [0] is the translated start and [1]
		 *         the end of location a, and [2] and [3] are the translated
		 *         start and end of location b.
		 */
		public static int[] modulateCircularLocationPair(Location a,
				Location b, int seqLength) {
			// Dummy case for non-circular locations.
			if (seqLength == 0)
				return new int[] { a.getMin(), a.getMax(), b.getMin(),
						b.getMax() };
			// Modulate our start/end to shortest possible equivalent region
			int[] aParts = modulateCircularLocation(a.getMin(), a.getMax(),
					seqLength);
			int aStart = aParts[0];
			int aEnd = aParts[1];
			// Modulate their start/end to shortest possible equivalent region
			int[] bParts = modulateCircularLocation(b.getMin(), b.getMax(),
					seqLength);
			int bStart = bParts[0];
			int bEnd = bParts[1];
			// If we wrap and the point we are checking for is before our start,
			// increment it by circularLength length
			if (aEnd > seqLength && bStart < aStart) {
				bStart += seqLength;
				bEnd += seqLength;
			}
			return new int[] { aStart, aEnd, bStart, bEnd };
		}

		/**
		 * Takes a point on a circular location and moves it left until it falls
		 * at the earliest possible point that represents the same base.
		 * 
		 * @param index
		 *            the point on the location to shift
		 * @param seqLength
		 *            the size of the circular location
		 * @return the shifted point
		 */
		public static int modulateCircularIndex(int index, int seqLength) {
			// Dummy case
			if (seqLength == 0)
				return index;
			// Modulate
			while (index > seqLength)
				index -= seqLength;
			return index;
		}

		/**
		 * Attempts to convert a plain Location into a RichLocation.
		 * 
		 * @param l
		 *            the location to convert
		 * @return the converted location
		 */
		public static RichLocation enrich(Location l) {
			// Dummy case where location is already enriched
			if (l instanceof RichLocation) {
				return (RichLocation) l;
			}
			// Compound case
			else if (l instanceof MergeLocation || !l.isContiguous()) {
				List<Location> members = new ArrayList<Location>();
				for (Iterator<Location> i = l.blockIterator(); i.hasNext();) {
					Location member = i.next();
					members.add(enrich(member));
				}
				return RichLocation.Tools.construct(RichLocation.Tools
						.merge(members));
			}
			// Fuzzy single points
			else if (l instanceof FuzzyPointLocation) {
				FuzzyPointLocation f = (FuzzyPointLocation) l;
				Position pos = new SimplePosition(f.hasBoundedMin(), f
						.hasBoundedMax(), f.getMin(), f.getMax(),
						Position.IN_RANGE);
				return new SimpleRichLocation(pos, 0); // 0 for no rank
			}
			// Fuzzy ranges
			else if (l instanceof FuzzyLocation) {
				FuzzyLocation f = (FuzzyLocation) l;
				Position start = new SimplePosition(f.hasBoundedMin(), false, f
						.getMin());
				Position end = new SimplePosition(false, f.hasBoundedMax(), f
						.getMax());
				return new SimpleRichLocation(start, end, 0); // 0 for no rank
			}
			// Normal ranges
			else if (l instanceof RangeLocation) {
				RangeLocation r = (RangeLocation) l;
				Position start = new SimplePosition(false, false, r.getMin());
				Position end = new SimplePosition(false, false, r.getMax());
				return new SimpleRichLocation(start, end, 0); // 0 for no rank
			}
			// Normal points
			else if (l instanceof PointLocation) {
				PointLocation p = (PointLocation) l;
				Position pos = new SimplePosition(false, false, p.getMin());
				return new SimpleRichLocation(pos, 0); // 0 for no rank
			}
			// Empty locations
			else if (l.toString().equals("{}")) {
				return EMPTY_LOCATION;
			}
			// All other cases
			else {
				throw new IllegalArgumentException(
						"Unable to enrich locations of type " + l.getClass());
			}
		}
	}
}
