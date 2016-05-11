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
 * Created on DATE
 *
 */
package org.biojava.nbio.core.sequence;

/**
 * @author Scooter Willis
 *
 */
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.compound.RNACompoundSet;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.ProxySequenceReader;
import org.biojava.nbio.core.sequence.template.SequenceView;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine;
import org.biojava.nbio.core.sequence.views.ComplementSequenceView;
import org.biojava.nbio.core.sequence.views.ReversedSequenceView;

/**
 * RNASequence where RNACompoundSet are the allowed values
 * @author Scooter Willis <willishf at gmail dot com>
 */

public class RNASequence extends AbstractSequence<NucleotideCompound> {

		/**
		 * Create a RNA sequence from a String
		 * @param seqString
		 * @throws CompoundNotFoundException
		 */
	public RNASequence(String seqString) throws CompoundNotFoundException {
		super(seqString, RNACompoundSet.getRNACompoundSet());
	}

	/**
	 * Create a RNA aequence from a proxy reader
	 * @param proxyLoader
	 */
	public RNASequence(ProxySequenceReader<NucleotideCompound> proxyLoader) {
		super(proxyLoader, RNACompoundSet.getRNACompoundSet());
	}

	/**
	 * Create a RNA sequence from a string with a user defined RNA compound set
	 * @param seqString
	 * @param compoundSet
	 * @throws CompoundNotFoundException
	 */
	public RNASequence(String seqString, CompoundSet<NucleotideCompound> compoundSet) throws CompoundNotFoundException {
		super(seqString, compoundSet);
	}

	/**
	 * Create a RNA sequence from a proxy reader and user defined RNA compound set
	 * @param proxyLoader
	 * @param compoundSet
	 */
	public RNASequence(ProxySequenceReader<NucleotideCompound> proxyLoader,
			CompoundSet<NucleotideCompound> compoundSet) {
		super(proxyLoader, compoundSet);
	}

	/**
	 * Get reverse complement view of the sequence
	 * @return
	 */
	public SequenceView<NucleotideCompound> getReverseComplement() {
		return new ComplementSequenceView<>(getInverse());
	}

	/**
	 * Get the inverse view of the sequence. It is the reverse sequence from
	 * end to begin where use reverse could imply complement. Called getInverse()
	 * in the hopes of making less confusing.
	 * @return
	 */
	@Override
public SequenceView<NucleotideCompound> getInverse() {
		return new ReversedSequenceView<>(this);
	}

	/**
	 * Get the complement view of the RNA sequence
	 * @return
	 */
	public SequenceView<NucleotideCompound> getComplement() {
		return new ComplementSequenceView<>(this);
	}

	/**
	 * Get the ProteinSequence from the RNA sequence
	 * @return
	 */
	public ProteinSequence getProteinSequence() {
		return getProteinSequence(TranscriptionEngine.getDefault());
	}

	/**
	 * Get the ProteinSequene from the RNA sequence with user defined
	 * transcription engine
	 *
	 * @param engine
	 * @return
	 */
	public ProteinSequence getProteinSequence(TranscriptionEngine engine) {
		return (ProteinSequence)engine.getRnaAminoAcidTranslator().createSequence(this);
	}

	public double getGC() {
		throw new UnsupportedOperationException("Not supported yet");
	}
}
