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

package org.biojava3.alignment;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.AlignedSequence.Step;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.template.ProfileView;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.AmbiguityRNACompoundSet;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.RNACompoundSet;
import org.biojava3.core.sequence.location.template.Location;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;


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
		return toString(width, null, getIDFormat(), true, true, true, true, true,false);
	}

	@Override
	public String toString(StringFormat format) {
		switch (format) {
		case ALN:
		case CLUSTALW:
		default:
			return toString(60, String.format("CLUSTAL W MSA from BioJava%n%n"), getIDFormat() + "   ", false, true,
					true, false, true,false);
		case FASTA:
			return toString(60, null, ">%s%n", false, false, false, false, false,false);
		case GCG:
		case MSF:
			String idFormat = getIDFormat();
			StringBuilder header = new StringBuilder();
			header.append(String.format("MSA from BioJava%n%n MSF: %d  Type: %s  Check: %d ..%n%n", getLength(),
					getGCGType(), getGCGChecksum()));
			for (AlignedSequence<S, C> as : list) {
				header.append(String.format(" Name: " + idFormat + " Len: %d  Check: %4d  Weight: %.1f%n",
						as.getAccession(), getLength(), getGCGChecksum(as), 1.0f)); // TODO show weights in MSF header
			}
			header.append(String.format("%n//%n%n"));
			// TODO? convert gap characters to '.'
			return toString(50, header.toString(), idFormat, false, false, true, false, false,false);
		case PDBWEB:            
			String pdbFormat = "%s";
			return toString(60,null,pdbFormat, true,true,true,false,true,true );
		}
	}

	// method from Object

	@Override
	public String toString() {
		return toString(getLength(), null, null, false, false, false, false, false,false);
	}

	// method for Iterable

	@Override
	public Iterator<AlignedSequence<S, C>> iterator() {
		return list.iterator();
	}

	// helper methods

	// calculates GCG checksum for entire Profile
	private int getGCGChecksum() {
		int check = 0;
		for (AlignedSequence<S, C> as : list) {
			check += getGCGChecksum(as);
		}
		return check % 10000;
	}

	// calculates GCG checksum for a given Sequence
	private int getGCGChecksum(AlignedSequence<S, C> sequence) {
		String s = sequence.toString().toUpperCase();
		int count = 0, check = 0;
		for (int i = 0; i < s.length(); i++) {
			count++;
			check += count * s.charAt(i);
			if (count == 57) {
				count = 0;
			}
		}
		return check % 10000;
	}

	// determines GCG type
	private String getGCGType() {
		CompoundSet<C> cs = getCompoundSet();
		return (cs == DNACompoundSet.getDNACompoundSet() || cs == AmbiguityDNACompoundSet.getDNACompoundSet()) ? "D" :
			(cs == RNACompoundSet.getRNACompoundSet() || cs == AmbiguityRNACompoundSet.getRNACompoundSet()) ? "R" :
				"P";
	}

	// creates format String for accession IDs
	private String getIDFormat() {
		int length = 0;
		for (AlignedSequence<S, C> as : list) {
			length = Math.max(length, (as.getAccession() == null) ? 0 : as.getAccession().toString().length());
		}
		return (length == 0) ? null : "%-" + (length + 1) + "s";
	}

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

						printSequenceAlignmentWeb(s,list, counter, idFormat,  seqIndexPre,  seqIndexFormatPre,  seqIndexPost,  seqIndexFormatPost,  start,  end);
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

					if ( aligConservation && list.size() == 2){
						if ( counter == 1)
							printConservation(s,idFormat,seqIndexPad,seqIndexPre, start,end,list, webDisplay);
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

		if ( webDisplay &&  aligConservation && list.size() == 2){
			s.append("</pre></div>");
			s.append("			<div class=\"subText\">");
			s.append("			<b>Legend:</b>");
			s.append("			<span class=\"m\">Black</span> - identical residues |"); 
			s.append("			<span class=\"sm\">Pink</span> - similar residues | ");
			s.append("			<span class=\"qg\">Blue</span> - sequence mismatch |");
			s.append("			<span class=\"dm\">Brown</span> - insertion/deletion |");			       
			s.append("		</div>");
			s.append(String.format("%n"));

		}
		return s.toString();
	}

	private void printSequenceAlignmentWeb(StringBuilder s,
			List<AlignedSequence<S, C>> list2, int counter, String idFormat, boolean seqIndexPre, String seqIndexFormatPre, boolean seqIndexPost, String seqIndexFormatPost, int start, int end) {

		AlignedSequence<S,C> as1 = list.get(0);
		AlignedSequence<S,C> as2 = list.get(1);

		
		AlignedSequence<S,C> as = list.get(counter - 1);

		if (idFormat != null) {
			s.append(String.format(idFormat, as.getAccession()));
		}
		if (seqIndexPre) {
			s.append(String.format(seqIndexFormatPre, as.getSequenceIndexAt(start)));
		}

        String mySeq = as.getSubSequence(start, end).getSequenceAsString();

        String s1 = as1.getSubSequence(start, end).getSequenceAsString();
        String s2 = as2.getSubSequence(start, end).getSequenceAsString();
		
		for ( int i =0 ; i < s1.length(); i++){
			if (  i >= s2.length())
				break;
			if ( i >= mySeq.length())
				break;
			
			char c1 = s1.charAt(i);
			char c2 = s2.charAt(i);
			char c = mySeq.charAt(i);
			if ( c1 == c2) 
				s.append("<span class=\"m\">"+c+"</span>");                    			
			else if (  isSimilar(c1,c2))  
				s.append("<span class=\"sm\">"+c+"</span>");
			else  if ( c1 == '-' || c2 == '-')				
				s.append("<span class=\"dm\">"+c+"</span>");
			else 
				s.append("<span class=\"qg\">"+c+"</span>");
			
		}
		
		

		if (seqIndexPost) {
			s.append(String.format(seqIndexFormatPost, as.getSequenceIndexAt(end)));
		}
		s.append(String.format("%n"));

	}

	private void printConservation(StringBuilder s, String idFormat,
			int seqIndexPad, boolean seqIndexPre, int start, int end, List<AlignedSequence<S, C>> list2, boolean webDisplay) {
		if ( list.size() == 2) {

			AlignedSequence<S,C> as1 = list.get(0);
			AlignedSequence<S,C> as2 = list.get(1);



			if (idFormat != null) {
				AccessionID ac1 = as1.getAccession();
				String id1 = "null";
				if ( ac1 != null)
					id1 = ac1.getID();

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

				if ( webDisplay) {

					if ( c1 == c2) 
						s.append("<span class=\"m\">|</span>");                    			
					else if (  isSimilar(c1,c2))  
						s.append("<span class=\"sm\">.</span>");
					else  if ( c1 == '-' || c2 == '-')
						s.append("<span class=\"dm\"> </span>");					
					else 
						s.append("<span class=\"qg\">&nbsp;</span>");
				} else {
					if ( c1 == c2) 
						s.append("|");                    			
					else if (  isSimilar(c1,c2))  
						s.append(".");					
					else 	
						s.append(" ");
				}

			}

			s.append(String.format("%n"));

		}

	}

	protected static final SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getBlosum65();

	private boolean isSimilar(char c1, char c2) {

		AminoAcidCompoundSet set = AminoAcidCompoundSet.getAminoAcidCompoundSet();

		AminoAcidCompound aa1 = set.getCompoundForString(""+c1);
		AminoAcidCompound aa2 = set.getCompoundForString(""+c2);

		short val = matrix.getValue(aa1,aa2);
		if ( val > 0 )
			return true;

		return false;
	}

}
