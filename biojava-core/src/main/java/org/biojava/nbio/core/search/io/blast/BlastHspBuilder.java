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
package org.biojava.nbio.core.search.io.blast;

import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.biojava.nbio.core.alignment.SimpleAlignedSequence;
import org.biojava.nbio.core.alignment.SimpleSequencePair;
import org.biojava.nbio.core.alignment.template.AlignedSequence.Step;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * Designed by Paolo Pavan.
 * You may want to find my contacts on Github and LinkedIn for code info
 * or discuss major changes.
 * https://github.com/paolopavan
 *
 * @author Paolo Pavan
 */

public class BlastHspBuilder {
	//private static final org.slf4j.Logger logger = LoggerFactory.getLogger(BlastHspBuilder.class);

	private int hspNum;
	private double hspBitScore;
	private int hspScore;
	private double hspEvalue;
	private int hspQueryFrom;
	private int hspQueryTo;
	private int hspHitFrom;
	private int hspHitTo;
	private int hspQueryFrame;
	private int hspHitFrame;
	private int hspIdentity;
	private int hspPositive;
	private int hspGaps;
	private int hspAlignLen;
	private String hspQseq;
	private String hspHseq;
	private String hspIdentityString;
	private Double percentageIdentity;
	private Integer mismatchCount;

	public BlastHspBuilder() {
	}

	public BlastHspBuilder setHspNum(int hspNum) {
		this.hspNum = hspNum;
		return this;
	}

	public BlastHspBuilder setHspBitScore(double hspBitScore) {
		this.hspBitScore = hspBitScore;
		return this;
	}

	public BlastHspBuilder setHspScore(int hspScore) {
		this.hspScore = hspScore;
		return this;
	}

	public BlastHspBuilder setHspEvalue(double hspEvalue) {
		this.hspEvalue = hspEvalue;
		return this;
	}

	public BlastHspBuilder setHspQueryFrom(int hspQueryFrom) {
		this.hspQueryFrom = hspQueryFrom;
		return this;
	}

	public BlastHspBuilder setHspQueryTo(int hspQueryTo) {
		this.hspQueryTo = hspQueryTo;
		return this;
	}

	public BlastHspBuilder setHspHitFrom(int hspHitFrom) {
		this.hspHitFrom = hspHitFrom;
		return this;
	}

	public BlastHspBuilder setHspHitTo(int hspHitTo) {
		this.hspHitTo = hspHitTo;
		return this;
	}

	public BlastHspBuilder setHspQueryFrame(int hspQueryFrame) {
		this.hspQueryFrame = hspQueryFrame;
		return this;
	}

	public BlastHspBuilder setHspHitFrame(int hspHitFrame) {
		this.hspHitFrame = hspHitFrame;
		return this;
	}

	public BlastHspBuilder setHspIdentity(int hspIdentity) {
		this.hspIdentity = hspIdentity;
		return this;
	}

	public BlastHspBuilder setHspPositive(int hspPositive) {
		this.hspPositive = hspPositive;
		return this;
	}

	public BlastHspBuilder setHspGaps(int hspGaps) {
		this.hspGaps = hspGaps;
		return this;
	}

	public BlastHspBuilder setHspAlignLen(int hspAlignLen) {
		this.hspAlignLen = hspAlignLen;
		return this;
	}

	public BlastHspBuilder setHspQseq(String hspQseq) {
		this.hspQseq = hspQseq;
		return this;
	}

	public BlastHspBuilder setHspHseq(String hspHseq) {
		this.hspHseq = hspHseq;
		return this;
	}

	public BlastHspBuilder setHspIdentityString(String hspIdentityString) {
		this.hspIdentityString = hspIdentityString;
		return this;
	}

	public BlastHspBuilder setPercentageIdentity(Double percentageIdentity) {
		this.percentageIdentity = percentageIdentity;
		return this;
	}

	public BlastHspBuilder setMismatchCount(Integer mismatchCount) {
		this.mismatchCount = mismatchCount;
		return this;
	}

	public <S extends Sequence<C>,C extends Compound>
	BlastHsp<S,C> createBlastHsp(Function<String,S> buildSeq) {
		SequencePair<S,C> aln = buildSequencePair(buildSeq, hspQseq, hspHseq, hspIdentityString );
		return new BlastHsp<>(hspNum, hspBitScore, hspScore, hspEvalue, hspQueryFrom, hspQueryTo, hspHitFrom, hspHitTo,
				hspQueryFrame, hspHitFrame, hspIdentity, hspPositive, hspGaps, hspAlignLen, aln, percentageIdentity, mismatchCount);
	}
	
	private static <S extends Sequence<C>,C extends Compound>
	SequencePair<S, C> buildSequencePair(Function<String, S> buildSeq,
			String query, String hit, String gaps) {
		S q = buildSeq.apply(query);
		S h = buildSeq.apply(hit);
		List<Step> steps = getAlignmentsSteps(gaps);
		return new SimpleSequencePair<>(
				new SimpleAlignedSequence<>(q, steps),
				new SimpleAlignedSequence<>(h, steps));
	}

	private static List<Step> getAlignmentsSteps(String gapped){
		return gapped.chars()
				.mapToObj((chr) -> (char)chr == '-' ? Step.GAP : Step.COMPOUND )
				.collect(Collectors.toList());
	}

//	private static SequencePair<?, ?> buildSequencePair(String query, String hit, String gaps) {
//		String concat = query+hit;
//		concat = concat.replace("-", "");
//		
//		try {
//			if (concat.matches("^[ACTG]+$")) {
//				DNASequence q = new DNASequence(query, DNACompoundSet.getDNACompoundSet());
//				DNASequence h = new DNASequence(hit, DNACompoundSet.getDNACompoundSet());
//				List<Step> steps = getAlignmentsSteps(gaps);
//				return new SimpleSequencePair<>(
//						new SimpleAlignedSequence<>(q, steps),
//						new SimpleAlignedSequence<>(h, steps));
//			} else if (concat.matches("^[ACUG]+$")) {
//				RNASequence q = new RNASequence(query, RNACompoundSet.getRNACompoundSet());
//				RNASequence h = new RNASequence(hit, RNACompoundSet.getRNACompoundSet());
//				List<Step> steps = getAlignmentsSteps(gaps);
//				return new SimpleSequencePair<>(
//						new SimpleAlignedSequence<>(q, steps),
//						new SimpleAlignedSequence<>(h, steps));
//			} else { 
//				ProteinSequence q = new ProteinSequence(query, AminoAcidCompoundSet.getAminoAcidCompoundSet());
//				ProteinSequence h = new ProteinSequence(hit, AminoAcidCompoundSet.getAminoAcidCompoundSet());
//				List<Step> steps = getAlignmentsSteps(gaps);
//				return new SimpleSequencePair<>(
//						new SimpleAlignedSequence<>(q, steps),
//						new SimpleAlignedSequence<>(h, steps));
//			}
//		} catch (CompoundNotFoundException ex) {
//			logger.error("Unexpected error, could not find compound when creating Sequence object from Hsp", ex);
//		}
//		return null;
//	}
//
//
//	private static Class<? extends Sequence<?>> guessSequenceType(String gappedSequenceString) {
//		String sequenceString = gappedSequenceString.replace("-", "");
//
//		if (sequenceString.matches("^[ACTG]+$"))
//			return DNASequence.class;
//		else if (sequenceString.matches("^[ACUG]+$"))
//			return RNASequence.class;
//		else
//			return ProteinSequence.class;
//	}
}
