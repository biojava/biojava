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
 * Created on June 14, 2010
 * Author: Mark Chapman
 */

package org.biojava.nbio.alignment;

import org.biojava.nbio.alignment.template.AlignedSequence;
import org.biojava.nbio.alignment.template.AlignedSequence.Step;
import org.biojava.nbio.alignment.template.Profile;
import org.biojava.nbio.alignment.template.ProfileView;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.io.util.IOUtils;
import org.biojava.nbio.core.sequence.location.template.Location;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.util.*;


/**
 * Implements a data structure for the results of sequence alignment.  Every {@link List} returned is unmodifiable.
 *
 * @author Mark Chapman
 * @param <S> each element of the alignment {@link Profile} is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public class SimpleProfile<S extends Sequence<C>, C extends Compound> implements Profile<S, C> {

	private List<AlignedSequence<S, C>> list;
	private List<S> originals;
	private int length;

	/**
	 * Creates a pair profile for the given already aligned sequences.
	 *
	 * @param query the first sequence of the pair
	 * @param target the second sequence of the pair
	 * @throws IllegalArgumentException if sequences differ in size
	 */
	protected SimpleProfile(AlignedSequence<S, C> query, AlignedSequence<S, C> target) {
		if (query.getLength() != target.getLength()) {
			throw new IllegalArgumentException("Aligned sequences differ in size");
		}
		list = new ArrayList<AlignedSequence<S, C>>();
		list.add(query);
		list.add(target);
		list = Collections.unmodifiableList(list);
		originals = new ArrayList<S>();
		originals.add((S) query.getOriginalSequence());
		originals.add((S) target.getOriginalSequence());
		originals = Collections.unmodifiableList(originals);
		length = query.getLength();
	}

	/**
	 * Creates a profile from a single sequence.
	 *
	 * @param sequence sequence to seed profile
	 */
	public SimpleProfile(S sequence) {
		List<Step> s = new ArrayList<Step>();
		for (int i = 0; i < sequence.getLength(); i++) {
			s.add(Step.COMPOUND);
		}
		list = new ArrayList<AlignedSequence<S, C>>();
		list.add(new SimpleAlignedSequence<S, C>(sequence, s));
		list = Collections.unmodifiableList(list);
		originals = new ArrayList<S>();
		originals.add(sequence);
		originals = Collections.unmodifiableList(originals);
		length = sequence.getLength();
	}

	/**
	 * Creates a pair profile for the given sequences.
	 *
	 * @param query the first sequence of the pair
	 * @param target the second sequence of the pair
	 * @param sx lists whether the query sequence aligns a {@link Compound} or gap at each index of the alignment
	 * @param xb number of {@link Compound}s skipped in the query sequence before the aligned region
	 * @param xa number of {@link Compound}s skipped in the query sequence after the aligned region
	 * @param sy lists whether the target sequence aligns a {@link Compound} or gap at each index of the alignment
	 * @param yb number of {@link Compound}s skipped in the target sequence before the aligned region
	 * @param ya number of {@link Compound}s skipped in the target sequence after the aligned region
	 * @throws IllegalArgumentException if alignments differ in size or given sequences do not fit in alignments
	 */
	protected SimpleProfile(S query, S target, List<Step> sx, int xb, int xa, List<Step> sy, int yb, int ya) {
		if (sx.size() != sy.size()) {
			throw new IllegalArgumentException("Alignments differ in size");
		}
		list = new ArrayList<AlignedSequence<S, C>>();
		list.add(new SimpleAlignedSequence<S, C>(query, sx, xb, xa));
		list.add(new SimpleAlignedSequence<S, C>(target, sy, yb, ya));
		list = Collections.unmodifiableList(list);
		originals = new ArrayList<S>();
		originals.add(query);
		originals.add(target);
		originals = Collections.unmodifiableList(originals);
		length = sx.size();
	}

	/**
	 * Creates a pair profile for the given profiles.
	 *
	 * @param query the first profile of the pair
	 * @param target the second profile of the pair
	 * @param sx lists whether the query profile aligns a {@link Compound} or gap at each index of the alignment
	 * @param sy lists whether the target profile aligns a {@link Compound} or gap at each index of the alignment
	 * @throws IllegalArgumentException if alignments differ in size or given profiles do not fit in alignments
	 */
	protected SimpleProfile(Profile<S, C> query, Profile<S, C> target, List<Step> sx, List<Step> sy) {
		if (sx.size() != sy.size()) {
			throw new IllegalArgumentException("Alignments differ in size");
		}
		list = new ArrayList<AlignedSequence<S, C>>();
		for (AlignedSequence<S, C> s : query) {
			list.add(new SimpleAlignedSequence<S, C>(s, sx));
		}
		for (AlignedSequence<S, C> s : target) {
			list.add(new SimpleAlignedSequence<S, C>(s, sy));
		}
		list = Collections.unmodifiableList(list);
		originals = new ArrayList<S>();
		originals.addAll(query.getOriginalSequences());
		originals.addAll(target.getOriginalSequences());
		originals = Collections.unmodifiableList(originals);
		length = sx.size();
	}

	 /**
     * Creates a profile for the already aligned sequences.
     * @param alignedSequences the already aligned sequences
     * @throws IllegalArgument if aligned sequences differ in length or
     * collection is empty.
     */
    public SimpleProfile(Collection<AlignedSequence<S,C>> alignedSequences) {
        list = new ArrayList<AlignedSequence<S,C>>();
        originals = new ArrayList<S>();
        
        Iterator<AlignedSequence<S,C>> itr = alignedSequences.iterator();
        if(!itr.hasNext()) {
            throw new IllegalArgumentException("alignedSequences must not be empty");
        }
        
        AlignedSequence<S, C> curAlignedSeq = itr.next();
        length = curAlignedSeq.getLength();
        list.add(curAlignedSeq);
        originals.add((S) curAlignedSeq.getOriginalSequence());
        
        while (itr.hasNext()) {
            curAlignedSeq = itr.next();
            if (curAlignedSeq.getLength() != length) {
                throw new IllegalArgumentException("Aligned sequences differ in size");
            }
            list.add(curAlignedSeq);
            originals.add((S) curAlignedSeq.getOriginalSequence());
        }
        list = Collections.unmodifiableList(list);
        originals = Collections.unmodifiableList(originals);
    }
	
	
	// methods for Profile

	@Override
	public AlignedSequence<S, C> getAlignedSequence(int listIndex) {
		return list.get(listIndex - 1);
	}

	@Override
	public AlignedSequence<S, C> getAlignedSequence(S sequence) {
		for (AlignedSequence<S, C> s : list) {
			if (s.equals(sequence) || s.getOriginalSequence().equals(sequence)) {
				return s;
			}
		}
		return null;
	}

	@Override
	public List<AlignedSequence<S, C>> getAlignedSequences() {
		return list;
	}

	@Override
	public List<AlignedSequence<S, C>> getAlignedSequences(int... listIndices) {
		List<AlignedSequence<S, C>> tempList = new ArrayList<AlignedSequence<S, C>>();
		for (int i : listIndices) {
			tempList.add(getAlignedSequence(i));
		}
		return Collections.unmodifiableList(tempList);
	}

	@Override
	public List<AlignedSequence<S, C>> getAlignedSequences(S... sequences) {
		List<AlignedSequence<S, C>> tempList = new ArrayList<AlignedSequence<S, C>>();
		for (S s : sequences) {
			tempList.add(getAlignedSequence(s));
		}
		return Collections.unmodifiableList(tempList);
	}

	@Override
	public C getCompoundAt(int listIndex, int alignmentIndex) {
		return getAlignedSequence(listIndex).getCompoundAt(alignmentIndex);
	}

	@Override
	public C getCompoundAt(S sequence, int alignmentIndex) {
		AlignedSequence<S, C> s = getAlignedSequence(sequence);
		return (s == null) ? null : s.getCompoundAt(alignmentIndex);
	}

	@Override
	public int[] getCompoundCountsAt(int alignmentIndex) {
		return getCompoundCountsAt(alignmentIndex, getCompoundSet().getAllCompounds());
	}

	@Override
	public int[] getCompoundCountsAt(int alignmentIndex, List<C> compounds) {
		int[] counts = new int[compounds.size()];
		C gap = getCompoundSet().getCompoundForString("-");
		int igap = compounds.indexOf(gap);
		for (C compound : getCompoundsAt(alignmentIndex)) {
			int i = compounds.indexOf(compound);
			if (i >= 0 && i != igap && !getCompoundSet().compoundsEquivalent(compound, gap)) {
				counts[i]++;
			}
		}
		return counts;
	}

	@Override
	public List<C> getCompoundsAt(int alignmentIndex) {
		// TODO handle circular alignments
		List<C> column = new ArrayList<C>();
		for (AlignedSequence<S, C> s : list) {
			column.add(s.getCompoundAt(alignmentIndex));
		}
		return Collections.unmodifiableList(column);
	}

	@Override
	public CompoundSet<C> getCompoundSet() {
		return list.get(0).getCompoundSet();
	}

	@Override
	public float[] getCompoundWeightsAt(int alignmentIndex) {
		return getCompoundWeightsAt(alignmentIndex, getCompoundSet().getAllCompounds());
	}

	@Override
	public float[] getCompoundWeightsAt(int alignmentIndex, List<C> compounds) {
		float[] weights = new float[compounds.size()];
		int[] counts = getCompoundCountsAt(alignmentIndex, compounds);
		float total = 0.0f;
		for (int i : counts) {
			total += i;
		}
		if (total > 0.0f) {
			for (int i = 0; i < weights.length; i++) {
				weights[i] = counts[i]/total;
			}
		}
		return weights;
	}

	@Override
	public int getIndexOf(C compound) {
		for (int i = 1; i <= length; i++) {
			if (getCompoundsAt(i).contains(compound)) {
				return i;
			}
		}
		return -1;
	}

	@Override
	public int[] getIndicesAt(int alignmentIndex) {
		int[] indices = new int[list.size()];
		for (int i = 0; i < indices.length; i++) {
			indices[i] = list.get(i).getSequenceIndexAt(alignmentIndex);
		}
		return indices;
	}

	@Override
	public int getLastIndexOf(C compound) {
		for (int i = length; i >= 1; i--) {
			if (getCompoundsAt(i).contains(compound)) {
				return i;
			}
		}
		return -1;
	}

	@Override
	public int getLength() {
		return length;
	}

	@Override
	public List<S> getOriginalSequences() {
		return originals;
	}

	@Override
	public int getSize() {
		int size = 0;
		for (AlignedSequence<S, C> s : list) {
			size += s.getOverlapCount();
		}
		return size;
	}

	@Override
	public ProfileView<S, C> getSubProfile(Location location) {
		// TODO ProfileView<S, C> getSubProfile(Location)
		return null;
	}

	@Override
	public boolean hasGap(int alignmentIndex) {
		C gap = getCompoundSet().getCompoundForString("-");
		for (C compound : getCompoundsAt(alignmentIndex)) {
			if (getCompoundSet().compoundsEquivalent(compound, gap)) {
				return true;
			}
		}
		return false;
	}

	@Override
	public boolean isCircular() {
		for (AlignedSequence<S, C> s : list) {
			if (s.isCircular()) {
				return true;
			}
		}
		return false;
	}

	@Override
	public String toString(int width) {
		return toString(width, null, IOUtils.getIDFormat(list), true, true, true, true, true, false);
	}

	@Override
	public String toString(StringFormat format) {
        switch (format) {
        case ALN:
        case CLUSTALW:
        default:
            return toString(60, String.format("CLUSTAL W MSA from BioJava%n%n"), IOUtils.getIDFormat(list) + "   ",
                    false, true, true, false, true, false);
        case FASTA:
            return toString(60, null, ">%s%n", false, false, false, false, false, false);
        case GCG:
        case MSF:
            return toString(50, IOUtils.getGCGHeader(list), IOUtils.getIDFormat(list), false, false, true, false,
                    false, false);
        case PDBWEB:
            return toString(60, null, "%10s", true, true, true, false, true, true);
        }
	}

	// method from Object

	@Override
	public String toString() {
		return toString(getLength(), null, null, false, false, false, false, false, false);
	}

	// method for Iterable

	@Override
	public Iterator<AlignedSequence<S, C>> iterator() {
		return list.iterator();
	}

	// helper methods

	// creates formatted String
	private String toString(int width, String header, String idFormat, boolean seqIndexPre, boolean seqIndexPost,
			boolean interlaced, boolean aligIndices, boolean aligConservation, boolean webDisplay) {

		// TODO handle circular alignments
		StringBuilder s = (header == null) ? new StringBuilder() : new StringBuilder(header);
		
		if ( webDisplay && list.size() == 2){
			s.append("<div><pre>");
		}

		width = Math.max(1, width);
		int seqIndexPad = (int) (Math.floor(Math.log10(getLength())) + 2);
		String seqIndexFormatPre = "%" + seqIndexPad + "d ", seqIndexFormatPost = "%" + seqIndexPad + "d";
		if (interlaced) {
			String aligIndFormat = "%-" + Math.max(1, width / 2) + "d %" + Math.max(1, width - (width / 2) - 1) +
			"d%n";
			for (int i = 0; i < getLength(); i += width) {
				int start = i + 1, end = Math.min(getLength(), i + width);
				if (i > 0) {
					s.append(String.format("%n"));
				}
				if (aligIndices) {
					if (end < i + width) {
						int line = end - start + 1;
						aligIndFormat = "%-" + Math.max(1, line / 2) + "d %" + Math.max(1, line - (line / 2) - 1) +
						"d%n";
					}
					if (idFormat != null) {
						s.append(String.format(idFormat, ""));
					}
					if (seqIndexPre) {
						s.append(String.format("%" + (seqIndexPad + 1) + "s", ""));
					}
					s.append(String.format(aligIndFormat, start, end));
				}

				int counter = 0;
				for (AlignedSequence<S, C> as : list) {
					counter++;

					if ( webDisplay && list.size() == 2){
						printSequenceAlignmentWeb(s, counter, idFormat, seqIndexPre, seqIndexFormatPre, seqIndexPost,
						        seqIndexFormatPost, start, end);
					} else {
						if (idFormat != null) {
							s.append(String.format(idFormat, as.getAccession()));
						}
						if (seqIndexPre) {
							s.append(String.format(seqIndexFormatPre, as.getSequenceIndexAt(start)));
						}

						s.append(as.getSubSequence(start, end).getSequenceAsString());

						if (seqIndexPost) {
							s.append(String.format(seqIndexFormatPost, as.getSequenceIndexAt(end)));
						}
						s.append(String.format("%n"));
					}

					if (aligConservation && list.size() == 2 && counter == 1) {
					    printConservation(s, idFormat, seqIndexPad, seqIndexPre, start, end, webDisplay);
					}
				}

			}
		} else {
			for (AlignedSequence<S, C> as : list) {
				if (idFormat != null) {
					s.append(String.format(idFormat, as.getAccession()));
				}
				for (int i = 0; i < getLength(); i += width) {
					int start = i + 1, end = Math.min(getLength(), i + width);
					if (seqIndexPre) {
						s.append(String.format(seqIndexFormatPre, as.getSequenceIndexAt(start)));
					}
					s.append(as.getSubSequence(start, end).getSequenceAsString());
					if (seqIndexPost) {
						s.append(String.format(seqIndexFormatPost, as.getSequenceIndexAt(end)));
					}
					s.append(String.format("%n"));
				}
			}
		}

		if (webDisplay && aligConservation && list.size() == 2) {
			s.append(IOUtils.getPDBLegend());
		}
		return s.toString();
	}

	private void printSequenceAlignmentWeb(StringBuilder s, int counter, String idFormat, boolean seqIndexPre,
	        String seqIndexFormatPre, boolean seqIndexPost, String seqIndexFormatPost, int start, int end) {
		AlignedSequence<S,C> as1 = list.get(0), as2 = list.get(1), as = list.get(counter - 1);

		if (idFormat != null) {
			s.append(String.format(idFormat, as.getAccession()));
		}
		if (seqIndexPre) {
			s.append(String.format(seqIndexFormatPre, as.getSequenceIndexAt(start)));
		}

        String mySeq = as.getSubSequence(start, end).getSequenceAsString();
        String s1 = as1.getSubSequence(start, end).getSequenceAsString();
        String s2 = as2.getSubSequence(start, end).getSequenceAsString();

		for (int i = 0; i < s1.length(); i++) {
			if (i >= s2.length() || i >= mySeq.length())
				break;

			char c1 = s1.charAt(i);
			char c2 = s2.charAt(i);
			char c = mySeq.charAt(i);
			s.append(IOUtils.getPDBCharacter(true, c1, c2, isSimilar(c1, c2), c));
		}

		if (seqIndexPost) {
			s.append(String.format(seqIndexFormatPost, as.getSequenceIndexAt(end)));
		}
		s.append(String.format("%n"));

	}

	private void printConservation(StringBuilder s, String idFormat, int seqIndexPad, boolean seqIndexPre, int start,
	        int end, boolean webDisplay) {
	    AlignedSequence<S,C> as1 = list.get(0), as2 = list.get(1);

	    if (idFormat != null) {
	        AccessionID ac1 = as1.getAccession();
	        String id1 = (ac1 == null) ? "null" : ac1.getID();
	        id1 = id1.replaceAll(".", " ");
	        s.append(String.format(idFormat, id1));
	    }

	    if (seqIndexPre) {
	        s.append(String.format("%" + (seqIndexPad + 1) + "s", ""));
	    }

	    String subseq1 = as1.getSubSequence(start, end).getSequenceAsString();
	    String subseq2 = as2.getSubSequence(start, end).getSequenceAsString();

	    for ( int ii =0 ; ii < subseq1.length() ; ii++){
	        if ( ii >= subseq2.length())
	            break;
	        char c1 = subseq1.charAt(ii);
	        char c2 = subseq2.charAt(ii);
	        s.append(IOUtils.getPDBConservation(webDisplay, c1, c2, isSimilar(c1, c2)));
	    }

	    s.append(String.format("%n"));
	}

	protected static final SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getBlosum65();

	private boolean isSimilar(char c1, char c2) {
		AminoAcidCompoundSet set = AminoAcidCompoundSet.getAminoAcidCompoundSet();

		AminoAcidCompound aa1 = set.getCompoundForString(""+c1);
		AminoAcidCompound aa2 = set.getCompoundForString(""+c2);

		short val = matrix.getValue(aa1,aa2);
		return val > 0;
	}

}
