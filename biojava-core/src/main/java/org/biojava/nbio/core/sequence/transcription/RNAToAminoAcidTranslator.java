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
package org.biojava.nbio.core.sequence.transcription;

import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava.nbio.core.sequence.template.AbstractCompoundTranslator;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.sequence.template.SequenceView;
import org.biojava.nbio.core.sequence.transcription.Table.Codon;
import org.biojava.nbio.core.sequence.views.WindowedSequence;

import java.util.*;

/**
 * Takes a {@link Sequence} of {@link NucleotideCompound} which should represent
 * an RNA sequence ({@link RNASequence} is good for this) and returns a list of
 * {@link Sequence} which hold {@link AminoAcidCompound}. The translator can
 * also trim stop codons as well as changing any valid start codon to an
 * initiating met.
 *
 * @author ayates
 */
public class RNAToAminoAcidTranslator extends
		AbstractCompoundTranslator<NucleotideCompound, AminoAcidCompound> {

	private final boolean trimStops;
	private final boolean initMetOnly;
	private final Map<Table.CaseInsensitiveTriplet, Codon> quickLookup;
	private final Map<AminoAcidCompound, List<Codon>> aminoAcidToCodon;
	// Cheeky lookup which uses a hashing value; key is to switch to using this
	// all the time
	private final Codon[] codonArray = new Codon[64000];
	private final AminoAcidCompound unknownAminoAcidCompound;
	private final AminoAcidCompound methionineAminoAcidCompound;
	private final boolean translateNCodons;

	// If true, then translation will stop at the first stop codon encountered
	// in the reading frame (the stop codon will be included as the last residue
	// in the resulting ProteinSequence, unless removed by #trimStops)
	private final boolean stopAtStopCodons;

	// If true, then translation will not start until the first start codon
	// encountered in the reading frame. The start codon will be included as the
	// first residue in the resulting ProteinSequence
	private final boolean waitForStartCodon;


	public RNAToAminoAcidTranslator(
			SequenceCreatorInterface<AminoAcidCompound> creator,
			CompoundSet<NucleotideCompound> nucleotides,
			CompoundSet<Codon> codons,
			CompoundSet<AminoAcidCompound> aminoAcids, Table table,
			boolean trimStops, boolean initMetOnly, boolean translateNCodons,
			boolean stopAtStopCodons, boolean waitForStartCodon) {

		super(creator, nucleotides, aminoAcids);
		this.trimStops = trimStops;
		this.initMetOnly = initMetOnly;
		this.translateNCodons = translateNCodons;

		quickLookup = new HashMap<Table.CaseInsensitiveTriplet, Codon>(codons
				.getAllCompounds().size());
		aminoAcidToCodon = new HashMap<AminoAcidCompound, List<Codon>>();

		List<Codon> codonList = table.getCodons(nucleotides, aminoAcids);
		for (Codon codon : codonList) {
			quickLookup.put(codon.getTriplet(), codon);
			codonArray[codon.getTriplet().intValue()] = codon;

			List<Codon> codonL = aminoAcidToCodon.get(codon.getAminoAcid());
			if (codonL == null) {
				codonL = new ArrayList<Codon>();
				aminoAcidToCodon.put(codon.getAminoAcid(), codonL);
			}
			codonL.add(codon);

		}
		unknownAminoAcidCompound = aminoAcids.getCompoundForString("X");
		methionineAminoAcidCompound = aminoAcids.getCompoundForString("M");
		this.stopAtStopCodons = stopAtStopCodons;
		this.waitForStartCodon = waitForStartCodon;
	}

	/**
	 * Performs the core conversion of RNA to Peptide. It does this by walking a
	 * windowed version of the given sequence. Any trailing DNA base pairs are
	 * ignored according to the specification of {@link WindowedSequence}.
	 */

	@Override
	public List<Sequence<AminoAcidCompound>> createSequences(
			Sequence<NucleotideCompound> originalSequence) {

		List<List<AminoAcidCompound>> workingList = new ArrayList<List<AminoAcidCompound>>();

		Iterable<SequenceView<NucleotideCompound>> iter = new WindowedSequence<NucleotideCompound>(
				originalSequence, 3);

		boolean first = true;

		// If not waiting for a start codon, start translating immediately
		boolean doTranslate = !waitForStartCodon;

		for (SequenceView<NucleotideCompound> element : iter) {
			AminoAcidCompound aminoAcid = null;

			int i = 1;
			Table.CaseInsensitiveTriplet triplet = new Table.CaseInsensitiveTriplet(
					element.getCompoundAt(i++), element.getCompoundAt(i++),
					element.getCompoundAt(i++));

			Codon target = null;

			target = quickLookup.get(triplet);

			// Check for a start
			if (!doTranslate && target.isStart()) {
				doTranslate = true;
			}

			if (doTranslate) {
				if (target != null)
					aminoAcid = target.getAminoAcid();
				if (aminoAcid == null && translateNCodons()) {
					aminoAcid = unknownAminoAcidCompound;
				} else {
					if (first && initMetOnly && target.isStart()) {
						aminoAcid = methionineAminoAcidCompound;
					}
				}

				addCompoundsToList(Arrays.asList(aminoAcid), workingList);
			}

			if (doTranslate && stopAtStopCodons && target.isStop()) {
				// Check if we need to stop, but dont stop until started!
				break;
			}

			first = false;
		}
		postProcessCompoundLists(workingList);

		return workingListToSequences(workingList);
	}

	/**
	 * Performs the trimming of stop codons and the conversion of a valid start
	 * amino acid to M
	 */
	@Override
	protected void postProcessCompoundLists(
			List<List<AminoAcidCompound>> compoundLists) {
		for (List<AminoAcidCompound> compounds : compoundLists) {
			if (trimStops) {
				trimStop(compounds);
			}
		}
	}

	/**
	 * Imperfect code. Checks the last amino acid to see if a codon could have
	 * translated a stop for it. Left in for the moment
	 */
	protected void trimStop(List<AminoAcidCompound> sequence) {
		AminoAcidCompound stop = sequence.get(sequence.size() - 1);
		boolean isStop = false;
		if (aminoAcidToCodon.containsKey(stop)) {
			for (Codon c : aminoAcidToCodon.get(stop)) {
				if (c.isStop()) {
					isStop = true;
					break;
				}
			}
		}

		if (isStop) {
			sequence.remove(sequence.size() - 1);
		}
	}

	/**
	 * Indicates if we want to force exact translation of compounds or not i.e.
	 * those with internal N RNA bases. This will cause a translation to an X
	 * amino acid
	 */
	public boolean translateNCodons() {
		return translateNCodons;
	}
}
