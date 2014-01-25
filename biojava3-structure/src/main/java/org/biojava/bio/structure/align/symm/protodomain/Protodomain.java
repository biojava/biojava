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
 * Created on 2012-11-20
 *
 */

package org.biojava.bio.structure.align.symm.protodomain;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomPositionMap;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.ResidueRange;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureIdentifier;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopDomain;

/**
 * A symmetry subunit for of some {@link StructureIdentifier enclosing structure}.
 * Responsible for defining the symmetry subunit of a structure, and for determining a unique name for each
 * Protodomain. This class is useful for requiring that a symmetry subunit created from an alignment remains within the
 * confines of its enclosing structure. A Protodomain created from an AFPChain is guaranteed to have this property.
 * To get the protodomain structure definition using CE-Symm, do:
 * <code>
 * 		AFPChain afpChain = ceSymm.align(ca1, ca2);
		int order = CeSymm.getSymmetryOrder(afpChain)
		String s = Protodomain.fromSymmetryAlignment(afpChain, ca1, order, atomCache).toString();
 * </code>
 * To get the protdomain structure string by using CE to align a structure against a protodomain of known
 * symmetry (and symmetry order), do:
 * <code>
 * 		AFPChain afpChain = ceSymm.align(ca1, ca2);
		String s = Protodomain.fromReferral(afpChain, ca2, order, atomCache).toString(); // note that we pass in ca2
 * </code>
 * 
 * @author dmyersturnbull
 * @since 3.08
 */
public class Protodomain implements StructureIdentifier {

	/**
	 * A helper class to create {@link ResidueRange ResidueRanges} without making a mistake.
	 * 
	 * @author dmyersturnbull
	 */
	private static class RangeBuilder {
		private Character currentChain;
		private ResidueNumber currentStart;
		private List<ResidueRange> list = new ArrayList<ResidueRange>();

		void addChain(char chainId) {
			if (currentChain != null) throw new IllegalStateException();
			currentChain = chainId;
		}

		void addEnd(ResidueNumber end, int length) {
			if (currentStart == null) throw new IllegalStateException();
			list.add(new ResidueRange(currentChain, currentStart, end, length));
			clearCurrent();
		}

		void addStart(ResidueNumber start) {
			if (currentChain == null) throw new IllegalStateException();
			if (currentStart != null) throw new IllegalStateException();
			currentStart = start;
		}

		void clearCurrent() {
			currentChain = null;
			currentStart = null;
		}

		List<ResidueRange> getList() {
			return list;
		}
	}

	/**
	 * If two {@link ResidueRange ResidueRanges} consecutive in the {@link #getResidueRanges() list of ranges} differ by no
	 * more than this number, they will be spliced together.
	 */
	private static final int APPROX_CONSECUTIVE = 4;

	private AtomCache cache;

	private int length;

	private final List<ResidueRange> list;

	private final String ranges;

	private final StructureIdentifier enclosingName;


	/**
	 * Client code should avoid calling this method if either fromSymmetryAlignment or fromReferral is more appropriate.
	 * The created Protodomain will include a gap in the alignment if and only if it is no greater than {@link #APPROX_CONSECUTIVE} residues.
	 * That is, the Protodomain will have two {@link ResidueRange ResidueRanges} for every gap that contains {@code APPROX_CONSECUTIVE} residues.
	 * Currently, {@code APPROX_CONSECUTIVE} is 4; a subsequent change will break the interface.
	 * 
	 * @param enclosing
	 *            The structure enclosing this Protodomain
	 * @param ranges
	 *            A list of Strings of the form {@code chain:} or {@code chain:start-end}
	 * @param afpChain
	 *            The AFPChain returned by the alignment. {@code afpChain.getName1()} and {@code afpChain.getName2()}
	 *            must be non-null and correct.
	 * @param ca
	 *            An array of C-alpha {@link Atom Atoms} from which to create the Protodomain
	 * @param keepStructure
	 *            If set to true, the structure will be recreated.
	 * @param numBlocks
	 *            The number of blocks in the alignment to include. Should either be 1 or 2.
	 * @param chainIndex
	 *            The chain in {@code afpChain} to create the Protodomain from.
	 * @return
	 * @throws ProtodomainCreationException
	 */
	public static Protodomain fromAfpChain(StructureIdentifier enclosing, AFPChain afpChain, Atom[] ca,
			int numBlocks, int chainIndex, AtomCache cache) throws ProtodomainCreationException {

		// int numAtomsInBlock1Alignment = afpChain.getOptLen()[0];
		// if (numAtomsInBlock1Alignment == 0) throw new ProtodomainCreationException("unknown", scopId);

		// this keeps a list of insertion-code-independent positions of the groups (where they are in the PDB file)
		final AtomPositionMap map = getAminoAcidPositions(cache, enclosing);

		// since CE-Symm has two blocks, we're going to have to do this in 2 parts:
		// from the start of block 1 to end end of block 1
		// PLUS from the start of block 2 to the end of block 2


		// we rely on the fact that SCOP won't give us a range that contains multiple chains
		// instead, any residues from another chain will be in a different range
		// however, CE-Symm CAN give us single blocks containing residues from different chains
		// we can deal with this by splitting each block into singleton sequential residues (e.g. A_7-7)
		// then we can splice these back together as we would with blocks

		final List<ResidueRange> totalRanges = new ArrayList<ResidueRange>();

		for (int block = 0; block < numBlocks; block++) {
			if (block + 1 > afpChain.getOptLen().length) {
				System.out.println("Warning: a block is missing in AFPChain from the alignment of " + enclosing
						+ " Assuming the alignment ends in block " + (block + 1) + ".");
				break;
			}
			final int numCaInEntireBlock = afpChain.getOptLen()[block];
			final int[] alignedPositions = afpChain.getOptAln()[block][chainIndex];
			// final int[] alignedPositions = afpChain.getOptAln()[block][chainIndex];
			for (int i = 0; i < numCaInEntireBlock; i++) {
				final int blockStart = alignedPositions[i];
				final int blockEnd = alignedPositions[i];
				final Group startGroup = ca[blockStart].getGroup();
				final Group endGroup = ca[blockEnd].getGroup();
				final ResidueNumber protodomainStartR = startGroup.getResidueNumber();
				final ResidueNumber protodomainEndR = endGroup.getResidueNumber();
				final int protodomainStart = map.getPosition(protodomainStartR);
				final int protodomainEnd = map.getPosition(protodomainEndR);
				final List<ResidueRange> blockRanges = getBlockRanges(protodomainStart, protodomainEnd, enclosing.getRanges(),
						protodomainStartR, protodomainEndR, startGroup, endGroup, map, enclosing);
				if (blockRanges != null) { // should this always be true?
					totalRanges.addAll(blockRanges);
				}
			}
		}

		final List<ResidueRange> splicedRanges;
		if (totalRanges.isEmpty()) {
			splicedRanges = totalRanges;
		} else {
			splicedRanges = spliceApproxConsecutive(map, totalRanges, APPROX_CONSECUTIVE);
		}

		final int length = ResidueRange.calcLength(splicedRanges);

		Protodomain protodomain = new Protodomain(enclosing, splicedRanges, length, cache);

		protodomain.length = length; // is the length of the AFPChain
		return protodomain;

	}

	/**
	 * Creates a Protodomain from an alignment between a Protodomain of known symmetry and some other structure.
	 * 
	 * @param afpChain
	 *            The <em>referred</em> ("result") structure must be in position 2, and the <em>referring</em> ("query") structure
	 *            in position 1.
	 * @param ca
	 *            An array of Atoms of the result
	 * @param cache
	 * @return
	 * @throws ProtodomainCreationException
	 */
	public static Protodomain fromReferral(StructureIdentifier enclosing, AFPChain afpChain, Atom[] ca, AtomCache cache)
			throws ProtodomainCreationException {
		if (afpChain.getOptLen().length != 1) throw new ProtodomainCreationException("unknown", afpChain.getName2(),
				"The AFPChain did not contain exactly 1 block.");
		try {
			return fromAfpChain(enclosing, afpChain, ca, 1, 1, cache);
		} catch (Exception e) { // IOException | StructureException | ArrayIndexOutOfBoundsException |
			// NullPointerException
			throw new ProtodomainCreationException("unknown", afpChain.getName2(), e);
		}
	}

	/**
	 * Creates a new Protodomain from a String. Ranges for this Protodomain will <em>not</em> be spliced;
	 * if this is desired, call {@link #spliceApproxConsecutive()} or {@link #spliceApproxConsecutive(int)}.
	 * 
	 * @param string
	 *            A string of the form: pdbId.chain_start-end,chain_start-end, ..., chain_start-end.
	 * @param enclosing
	 * @return
	 * @throws ProtodomainCreationException
	 */
	public static Protodomain fromString(String string, StructureIdentifier enclosing, AtomCache cache)
			throws ProtodomainCreationException, IllegalArgumentException {

		final AtomPositionMap map = Protodomain.getAminoAcidPositions(cache, enclosing);
		final List<ResidueRange> list = ResidueRange.parseMultiple(string.substring(string.indexOf('.') + 1), map);

		int length = ResidueRange.calcLength(list);

		Protodomain protodomain = new Protodomain(enclosing, list, length, cache);
		return protodomain;
	}

	/**
	 * Creates a Protodomain from an AFPChain returned by CE-Symm.
	 * 
	 * @param afpChain
	 * @param ca
	 *            An array of Atoms of the result
	 * @param cut
	 *            What the whole protodomain will be "cut" by. This should ordinarily be set to the symmetry order. If
	 *            it is set to 1, it the whole protodomain will be returned.
	 * @param cache
	 * @return
	 * @throws ProtodomainCreationException
	 */
	public static Protodomain fromSymmetryAlignment(StructureIdentifier enclosing, AFPChain afpChain, Atom[] ca, int order, AtomCache cache)
			throws ProtodomainCreationException {
		// note that order isn't required in the current implementation, but it could be
		if (afpChain.getBlockNum() != 2) throw new ProtodomainCreationException("unknown", afpChain.getName2(),
				"The AFPChain did not contain exactly 2 blocks.");
		return fromAfpChain(enclosing, afpChain, ca, 2, 0, cache);
	}

	/**
	 * Builds a new Protodomain from an array of protodomains of the same enclosing structure.
	 * @param map
	 * @param cache
	 * @param protodomains An array of protodomains, in order from first to last as defined in the enclosing structure
	 * @return
	 * TODO Test and make public
	 */
	private static Protodomain append(AtomPositionMap map, AtomCache cache, Protodomain... protodomains) {
		if (protodomains.length == 0) return null;
		if (protodomains.length == 1) return new Protodomain(protodomains[0]);
		StructureIdentifier scopId = protodomains[0].getEnclosingStructureIdentifier();
		List<ResidueRange> ranges = new ArrayList<ResidueRange>();
		int previousEnd = 0;
		for (Protodomain p : protodomains) {
			if (!p.getEnclosingStructureIdentifier().getIdentifier().equalsIgnoreCase(scopId.getIdentifier())) {
				throw new IllegalArgumentException("Protodomains do not belong to the same enclosing structure");
			}
			int start = map.getPosition(p.getResidueRanges().get(p.getResidueRanges().size()-1).getStart());
			if (start < previousEnd) {
				throw new IllegalArgumentException("Protodomains must be in order");
			}
			previousEnd = map.getPosition(p.getResidueRanges().get(p.getResidueRanges().size()-1).getEnd());
			ranges.addAll(p.getResidueRanges());
		}
		int length = ResidueRange.calcLength(ranges);
		Protodomain protodomain = new Protodomain(scopId, ranges, length, cache);
		return protodomain;
	}
	
	/**
	 * From a list of residue ranges {@code residueRanges} and a number of residues wanted {@code numResiduesWant} (
	 * <em>which includes any alignment gaps</em>), returns a new list of residue ranges that contains the first
	 * {@code numResiduesWant} residues from {@code residueRanges}.
	 * 
	 * @param residueRanges
	 * @param map
	 *            An AtomPositionMap for the enclosing structure
	 * @param order
	 *            The order of symmetry of this Protodomain
	 * @param index The index of the cut (sub-Protodomain): 0 will cause the first cut to be returned, 1 will cause the second cut to be returned, etc
	 * @return
	 */
	private static List<ResidueRange> calcSubstruct(List<ResidueRange> residueRanges, AtomPositionMap map, int order,
			int index) {

		if (order < 1) throw new IllegalArgumentException("Can't compute substructure if order is " + order);
		
		List<ResidueRange> part = new ArrayList<ResidueRange>(); // the parts we want to keep

		int numResiduesHave = 0;
		int numResiduesSkipped = 0;

		int numResiduesWant = ResidueRange.calcLength(residueRanges) / order;

		int numResiduesWantToSkip = numResiduesWant * index;
		if (index > 0) numResiduesWantToSkip++; // this is because we don't want overlap

		final NavigableMap<ResidueNumber, Integer> navMap = map.getNavMap();

		outer: for (ResidueRange rr : residueRanges) {

			if (numResiduesSkipped + rr.getLength() <= numResiduesWantToSkip) {
				numResiduesSkipped += rr.getLength();
				continue; // skip all of it
			} else if (numResiduesSkipped <= numResiduesWantToSkip) {
				// we only want part
				// we'll execute the block below
				// but before doing this, we'll redefine rr to remove the unwanted residues
				ResidueNumber redefStart = rr.getStart();
				for (int j = 0; j < numResiduesWantToSkip - numResiduesSkipped; j++) {
					redefStart = navMap.higherEntry(redefStart).getKey();
				}
				rr = new ResidueRange(rr.getChainId(), redefStart, rr.getEnd(), rr.getLength() - (numResiduesWantToSkip - numResiduesSkipped));
				numResiduesSkipped += rr.getLength();
			}
			
			// note that getLength() here DOES INCLUDE gaps (and it should)
			if (numResiduesHave + rr.getLength() <= numResiduesWant) {
				// no problem: we can fit the entire residue range in
				part.add(rr);

			} else if (numResiduesHave <= numResiduesWant) {
				// okay, so we can only fit some of rr in (and then we'll end)
				// we'll insert rr from rr.getStart() to rr.getStart() + numResiduesWant - numResiduesHave
				// unfortunately, this is harder because of insertion codes
				for (Map.Entry<ResidueNumber, Integer> entry : navMap.entrySet()) { // this is in
					// insertion-code-independent order
					if (entry.getKey().equals(rr.getStart())) {
						ResidueNumber endResidueNumber = entry.getKey(); // where we want to end this last residue range
						// okay, we got to the end of rr
						// now we need to go back by numResiduesWant - numResiduesTraversed
						for (int j = 0; j < numResiduesWant - numResiduesHave; j++) {
							endResidueNumber = navMap.higherEntry(endResidueNumber).getKey();
						}
						part.add(new ResidueRange(rr.getChainId(), rr.getStart(), endResidueNumber, numResiduesWant
								- numResiduesHave));
						break outer; // we filled up numResiduesTraversed, so we won't be able to add further residue
						// ranges
					}
				}
			}
			numResiduesHave += rr.getLength();
		}
		return part;
	}

	/**
	 * Really just calls {@link ResidueRange#getAminoAcidPositions(Atom[])}.
	 * 
	 * @param cache
	 * @param enclosing
	 * 			The enclosing structure of this Protodomain
	 * @return
	 * @throws ProtodomainCreationException
	 */
	private static AtomPositionMap getAminoAcidPositions(AtomCache cache, StructureIdentifier enclosing)
			throws ProtodomainCreationException {
		try {
			// We cannot use CA atoms only here because sometimes the C-alpha atom is missing
			// Our AtomPositionMap should use something more liberal (see the AtomPositionMap constructor)
			final Atom[] allAtoms = StructureTools.getAllAtomArray(cache.getStructure(enclosing.getPdbId()));
			// ok here?
			return new AtomPositionMap(allAtoms);
		} catch (IOException e) {
			throw new ProtodomainCreationException("unknown", enclosing.getIdentifier(), e,
					"Could not get a list of amino acid residue number positions.");
		} catch (StructureException e) {
			throw new ProtodomainCreationException("unknown", enclosing.getIdentifier(), e,
					"Could not get a list of amino acid residue number positions.");
		}
	}

	/**
	 * Finds the list of residue ranges that corresponds to the protodomain between {@code protodomainStartR} and
	 * {@code protodomainEndR} <em>within the {@link ScopDomain} identified by {@code domainRanges}</em>.
	 * 
	 * @param protodomainStart
	 *            The Protodomain's ATOM record start position
	 * @param protodomainEnd
	 *            The Protodomain's ATOM record end position
	 * @param domainRanges
	 *            A list of ranges returned by calling {@link ScopDomain#getRanges()} on this Protodomain's
	 *            {@link ScopDomain}.
	 * @param protodomainStartR
	 *            The Protodomain's ResidueNumber start position
	 * @param protodomainEndR
	 *            The Protodomain's ResidueNumber end position
	 * @param blockStartGroup
	 *            The {@link Group} at the start of a block in the the alignment
	 * @param blockEndGroup
	 *            The {@link Group} at the end of a block in the the alignment
	 * @param map
	 *         An AtomPositionMap for the enclosing structure
	 * @return A list of {@link ResidueRange ResidueRanges} within the overlap (intersection) of the relevant block of
	 *         the alignment corresponding to this Protodomain, and its enclosing structure's boundaries.
	 * @throws ProtodomainCreationException
	 */
	private static List<ResidueRange> getBlockRanges(int protodomainStart, int protodomainEnd,
			List<String> domainRanges, ResidueNumber protodomainStartR, ResidueNumber protodomainEndR,
			Group blockStartGroup, Group blockEndGroup, AtomPositionMap map, StructureIdentifier enclosing)
					throws ProtodomainCreationException {

//		if (domainRanges.isEmpty()) {
//			domainRanges.add(new ResidueRange(map.getFirst().getChainId(), map.getFirst(), map.getLast(), map.calcLength(map.getFirst(), map.getLast())).toString());
//		}
		
		/*
		 * Okay, things get annoying here, since we need to handle:
		 * 1. domain insertions, AND
		 * 2. heteroatom Groups inside this domain that are past the last amino acid
		 * If we simply take the protodomain's start and its end (in this block), we could include:
		 * 1. domains inserted between, AND
		 * 2. residues in other domains between the last amino acid in the correct domain and the last heteroatom (especially a water molecule) in the correct domain
		 * These are bad. So: we use SCOP's list of ranges for our domain, making sure that our protodomain does not extend outside
		 * its boundaries.
		 * To handle #2, we -could- get the last amino acid (for the end residue; the first aa for the start), but this is problematic for modified amino acids. We don't do this.
		 * Some heteroatoms also have alpha carbons.
		 * But that's not all. What about multi-chain domains? We need to handle those.
		 * Then things can get even MORE cumbersome thanks to insertion codes. d1qdmc1
		 * has an insertion code in one of its SCOP ranges: C:1S-104S (Residue numbers in ATOMs go: 247 1S 2S ... 103S 104S 248)
		 */

		if (protodomainStart > protodomainEnd) return null; // it's okay if our protodomain doesn't exist in this block

		RangeBuilder rangeBuilder = new RangeBuilder();
		for (String domainRange : domainRanges) {

			// a range is of the form: A:05-117 (chain:start-end) OR just A:
			ResidueRange rr = ResidueRange.parse(domainRange, map);
			ResidueNumber domainStartR = rr.getStart();
			ResidueNumber domainEndR = rr.getEnd();
			int domainStart, domainEnd;

			// these are the insertion-code-independent positions
			Integer domainStartO = map.getPosition(domainStartR);
			if (domainStartO == null) throw new ProtodomainCreationException("unknown", enclosing.getIdentifier(),
					"Couldn't find the start of the SCOP domain in the PDB file, which was supposed to be "
							+ domainStartR.printFull());
			domainStart = domainStartO;
			Integer domainEndO = map.getPosition(domainEndR);
			if (domainEndO == null) throw new ProtodomainCreationException("unknown", enclosing.getIdentifier(),
					"Couldn't find the end of the SCOP domain in the PDB file, which was supposed to be "
							+ domainEndR.printFull());
			domainEnd = domainEndO;

			// in this case, there's something wrong with the SCOP definition
			if (domainStart > domainEnd) return null;

			final char chain = domainRange.charAt(0);
			rangeBuilder.addChain(chain);

			// the domain part starts and ends before the start of the protodomain
			if (domainEnd < protodomainStart) {
				rangeBuilder.clearCurrent();
				continue; // well, obviously we're not using that part of the domain
				// the domain part starts before the protodomain starts but
			} else if (domainStart <= protodomainStart) {
				// ends after the protodomain starts
				rangeBuilder.addStart(protodomainStartR); // protodomain start
				// the domain part ends before the protodomain ends
				if (domainEnd <= protodomainEnd) {
					// domain end
					rangeBuilder.addEnd(domainEndR, map.calcLength(domainEnd, protodomainStart, chain));
					// the domain part ends after the protodomain ends
				} else {
					// protodomain end
					rangeBuilder.addEnd(protodomainEndR, map.calcLength(protodomainEnd, protodomainStart, chain));
				}
				// domain part ends before the protodomain starts
			} else if (domainStart > protodomainEnd) {
				rangeBuilder.clearCurrent();
				continue; // if we knew the domain parts are ordered, we could break
				// the domain part starts after the protodomain starts but ends
			} else if (domainStart > protodomainStart) {
				// after it starts
				rangeBuilder.addStart(domainStartR); // domain start
				if (domainEnd <= protodomainEnd) {
					rangeBuilder.addEnd(domainEndR, map.calcLength(domainEnd, domainStart, chain)); // domain
					// end
				} else {
					rangeBuilder.addEnd(protodomainEndR, map.calcLength(protodomainEnd, domainStart, chain)); // protodomain
					// end
				}
			} else { // this can't happen
				assert(false); // might as well catch a bug immediately
			}
		}

		return rangeBuilder.getList();
	}

	/**
	 * @see #spliceApproxConsecutive(int)
	 */
	private static List<ResidueRange> spliceApproxConsecutive(AtomPositionMap map, List<ResidueRange> old,
			int approxConsecutive) {
		List<ResidueRange> spliced = new ArrayList<ResidueRange>(old.size());

		ResidueRange growing = old.get(0);

		for (int i = 1; i < old.size(); i++) {

			final ResidueRange next = old.get(i);

			final int dist = map.calcLengthDirectional(growing.getEnd(), next.getStart());

			// if they're of different chains, CANNOT splice safely, so don't
			if (dist > approxConsecutive || !growing.getChainId().equals(next.getChainId())) {
				spliced.add(growing);
				growing = next;
			} else {
				growing = new ResidueRange(growing.getChainId(), growing.getStart(), next.getEnd(), growing.getLength()
						+ next.getLength() + dist);
			}

		}

		spliced.add(growing);

		return spliced;
	}

	/**
	 * @param pdbId
	 * @param enclosing
	 * @param list
	 *            A list of {@link ResidueRange ResidueRanges} that define this Protodomain.
	 * @param length
	 *            The number of amino acids contained in this Protodomain, <em>including any alignment gaps</em>. If
	 *            each ResidueRange in {@code list} has a non-null {@link ResidueRange#getLength() length}, the sum of
	 *            these lengths should be equal to this argument.
	 * @param cache
	 */
	public Protodomain(StructureIdentifier enclosing, List<ResidueRange> list, int length, AtomCache cache) {
		this.enclosingName = enclosing;
		this.list = list;
		this.cache = cache;
		this.length = length;
		StringBuilder sb = new StringBuilder();
		Iterator<ResidueRange> iter = list.iterator();
		while (iter.hasNext()) {
			ResidueRange t = iter.next();
			sb.append(t);
			if (iter.hasNext()) sb.append(",");
		}
		ranges = sb.toString();
	}

	public Protodomain(Protodomain protodomain) {
		this.cache = protodomain.cache;
		this.enclosingName = protodomain.enclosingName;
		this.length = protodomain.length;
		this.list = protodomain.list;
		this.ranges = protodomain.ranges;
	}

	/**
	 * Creates a new Protodomain corresponding to this Protodomain cut by {@code order}. For example, if this
	 * Protodomain contains 6 symmetry subunits, calling this.createSubgroup(3) will return the Protodomain
	 * corresponding to the first 2 symmetry subunits.
	 * 
	 * This method always returns the first cut (index 1) available. To get another cut, see {@link #createSubstruct(int, int)}.
	 * 
	 * @param order The order of symmetry of this Protodomain
	 * @return
	 * @throws ProtodomainCreationException
	 */
	public Protodomain createSubstruct(int order) throws ProtodomainCreationException {
		return createSubstruct(order, 0);
	}

	/**
	 * Creates a new Protodomain corresponding to this Protodomain cut by {@code order}. For example, if this
	 * Protodomain contains 6 symmetry subunits, calling this.createSubgroup(3) will return the Protodomain
	 * corresponding to the first 2 symmetry subunits.
	 * 
	 * @param order The order of symmetry of this Protodomain
	 * @param index The index of the cut (sub-Protodomain): 0 will cause the first cut to be returned, 1 will cause the second cut to be returned, etc
	 * @return
	 * @throws ProtodomainCreationException
	 * TODO Test and make public
	 */
	private Protodomain createSubstruct(int order, int index) throws ProtodomainCreationException {
		AtomPositionMap map = Protodomain.getAminoAcidPositions(cache, enclosingName);
		List<ResidueRange> rrs = calcSubstruct(list, map, order, index);
		return new Protodomain(enclosingName, rrs, ResidueRange.calcLength(rrs), cache);
	}


	public AtomCache getCache() {
		return cache;
	}

	/**
	 * @return The number of amino acids in this Protodomain, including any alignment gaps.
	 */
	public int getLength() {
		return length;
	}

	public String getPdbId() {
		return enclosingName.getPdbId();
	}

	/**
	 * @return The list of {@link ResidueRange ResidueRanges} in this Protodomain.
	 */
	public List<ResidueRange> getResidueRanges() {
		return list;
	}

	/**
	 * @return The list of ranges in this Protodomain.
	 */
	public List<String> getRanges() {
		return ResidueRange.toStrings(list);
	}

	/**
	 * @return A String describing the list of {@link ResidueRange ResidueRanges} of this Protodomain of the format
	 *         chain_start-end,chain.start-end, ... For example, <code>A_5-100,110-144</code>.
	 */
	public String getRangeString() {
		return ranges;
	}

	/**
	 * Returns this Protodomain's parent {@link StructureIdentifier}.
	 */
	public StructureIdentifier getEnclosingStructureIdentifier() {
		return enclosingName;
	}

	/**
	 * Sets the {@link AtomCache} used in this Protodomain. This would only need to be called if
	 * {@link #buildStructure()} will be called, and a different AtomCache is needed.
	 * 
	 * @param cache
	 */
	public void setCache(AtomCache cache) {
		this.cache = cache;
	}

	/**
	 * Splices residues that are approximately consecutive within a distance of {@link #APPROX_CONSECUTIVE}. This is the
	 * same value that is used internally for Protodomains created from AFPChains. Note that Protodomains created from
	 * Strings are never spliced during creation; manually call this method if that is desired.
	 * 
	 * @see #spliceApproxConsecutive(int)
	 */
	public Protodomain spliceApproxConsecutive() {
		return spliceApproxConsecutive(APPROX_CONSECUTIVE);
	}

	/**
	 * Splices residues that are approximately consecutive, within some distance. That is, from a list of
	 * {@link ResidueRange ResidueRanges}, returns a new list of ResidueRanges containing every range from {@code old},
	 * but with consecutive ranges starting and ending within {@code approxConsecutive} amino acids spliced together.
	 * Note that ranges from different chains are assigned infinite distance, so they cannot be spliced together. For
	 * example: 1dkl.A_1-100,A_102-200 would be spliced into 1dkl.A_1-200 for {@code approxConsecutive=2} or higher.
	 * 
	 * @param old
	 * @param posns
	 * @param sorted
	 * @param approxConsecutive
	 * @return
	 * @see {@link #APPROX_CONSECUTIVE}.
	 */
	public Protodomain spliceApproxConsecutive(int approxConsecutive) {
		try {
			List<ResidueRange> ranges = spliceApproxConsecutive(getAminoAcidPositions(cache, enclosingName), list,
					approxConsecutive);
			return new Protodomain(enclosingName, ranges, 1, cache);
		} catch (ProtodomainCreationException e) { // this really shouldn't happen, unless the AtomCache has changed
			// since this Protodomain was created
			throw new IllegalArgumentException(e);
		}
	}

	/**
	 * @return A String describing the list of {@link ResidueRange ResidueRanges} of this Protodomain of the format
	 *         pdbId.chain_start-end,chain.start-end, ... For example, <code>1w0p.A_5-100,110-144</code>.
	 */
	@Override
	public String toString() {
		return getIdentifier();
	}


	@Override
	public String getIdentifier() {
		if (ranges.isEmpty()) return enclosingName.getPdbId();
		return enclosingName.getPdbId() + "." + ranges;
	}


	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((getIdentifier() == null) ? 0 : getIdentifier().hashCode());
		return result;
	}


	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		Protodomain other = (Protodomain) obj;
		if (getIdentifier() == null) {
			if (other.getIdentifier() != null) return false;
		} else if (!getIdentifier().equals(other.getIdentifier())) return false;
		return true;
	}

}
