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
package org.biojava.nbio.core.sequence.transcription;

import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava.nbio.core.sequence.template.AbstractCompoundTranslator;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.ProxySequenceReader;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.sequence.views.RnaSequenceView;

import java.util.ArrayList;
import java.util.List;

/**
 * Performs the first stage of transcription by going from DNA to RNA. This
 * class will first delegate to {@link Frame} in order to be in the correctly
 * specified translation frame and then translates T to U. The other
 * translation carried out is to convert an equivalent compound in DNA to RNA
 * i.e. for the base A in DNA fetching the equivalent A base in the RNA
 * {@link CompoundSet}.
 *
 * @author ayates
 */
public class DNAToRNATranslator extends AbstractCompoundTranslator<NucleotideCompound, NucleotideCompound> {

		private final boolean shortCutTranslation;

	public DNAToRNATranslator(SequenceCreatorInterface<NucleotideCompound> rnaCreator,
			CompoundSet<NucleotideCompound> dna, CompoundSet<NucleotideCompound> rna,
			boolean shortCutTranslation) {
		super(rnaCreator, dna, rna);
		this.shortCutTranslation = shortCutTranslation;
		defaultMappings();
		thyamineToUracil();
	}

		/**
		 * Overloaded local version which delegates to an optional translator
		 * when told to (specified during construction).
		 *
		 * @param originalSequence The DNA sequence to translate
		 * @return The translated single sequence
		 */
		@Override
		public List<Sequence<NucleotideCompound>> createSequences(Sequence<NucleotideCompound> originalSequence) {
				if(shortCutTranslation) {
						List<Sequence<NucleotideCompound>> result = new ArrayList<Sequence<NucleotideCompound>>(1);
						result.add(wrapToRna(originalSequence));
						return result;
				}
				else {
						return super.createSequences(originalSequence);
				}
		}

		/**
		 * Takes in the given DNA Sequence and returns an instance of RNASequence
		 * which is using {@link RnaSequenceView} as a
		 * {@link ProxySequenceReader}.
		 */
		protected RNASequence wrapToRna(Sequence<NucleotideCompound> dna) {
				ProxySequenceReader<NucleotideCompound> rnaView = new RnaSequenceView(dna);
				return new RNASequence(rnaView);
		}

	private void defaultMappings() {
		NucleotideCompound thymine = getFromCompoundSet().getCompoundForString("T");
		for(NucleotideCompound dnaBase: getFromCompoundSet().getAllCompounds()) {
			if(dnaBase.equalsIgnoreCase(thymine)) {
				continue;
			}
			NucleotideCompound rnaBase = getToCompoundSet().getCompoundForString(
					dnaBase.toString());
			addCompounds(dnaBase, rnaBase);
		}

	}

	private void thyamineToUracil() {
		addCompounds(getFromCompoundSet().getCompoundForString("T"),
				getToCompoundSet().getCompoundForString("U"));
		addCompounds(getFromCompoundSet().getCompoundForString("t"),
				getToCompoundSet().getCompoundForString("u"));
	}

	public Sequence<NucleotideCompound> createSequence(Sequence<NucleotideCompound> originalSequence, Frame frame) {
		Sequence<NucleotideCompound> wrapped = frame.wrap(originalSequence);
		return super.createSequence(wrapped);
	}

	@Override
	public Sequence<NucleotideCompound> createSequence(Sequence<NucleotideCompound> originalSequence) {
		return createSequence(originalSequence, Frame.getDefaultFrame());
	}

	@Override
	protected void postProcessCompoundLists(
			List<List<NucleotideCompound>> compoundLists) {
		//No post processing needed
	}
}
