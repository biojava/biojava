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
package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.ProxySequenceReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Locale;

/**
 * A sequence creator which preserves the case of its input string in
 * the user collection of the returned ProteinSequence.
 *
 * <p>The user collection will be the same length as the resulting ProteinSequence.
 * Each object can be cast to a Boolean. If true, the corresponding position in
 * the input file was uppercase.
 *
 * <h3>Example</h3>
 * <code><pre>CasePreservingProteinSequenceCreator creator =
 *    new CasePreservingProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet());
 *AbstractSequence<AminoAcidCompound> seq = creator.getSequence("aaAA",0);
 *System.out.println(seq.getSequenceAsString()); //"AAAA"
 *System.out.println(seq.getUserCollection()); //"[false, false, true, true]"
 *</code></pre>
 */
public class CasePreservingProteinSequenceCreator extends ProteinSequenceCreator {

	private final static Logger logger = LoggerFactory.getLogger(CasePreservingProteinSequenceCreator.class);

	public CasePreservingProteinSequenceCreator(
			CompoundSet<AminoAcidCompound> compoundSet) {
		super(compoundSet);
	}

	/**
	 *
	 * @see org.biojava.nbio.core.sequence.io.ProteinSequenceCreator#getSequence(org.biojava.nbio.core.sequence.template.ProxySequenceReader, long)
	 */
	@Override
	public AbstractSequence<AminoAcidCompound> getSequence(
			ProxySequenceReader<AminoAcidCompound> proxyLoader, long index) {
		AbstractSequence<AminoAcidCompound> seq = super.getSequence(proxyLoader, index);
		seq.setUserCollection(getStringCase(proxyLoader.getSequenceAsString()));
		return seq;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.core.sequence.io.ProteinSequenceCreator#getSequence(java.lang.String, long)
	 */
	@Override
	public AbstractSequence<AminoAcidCompound> getSequence(String sequence,
			long index) throws CompoundNotFoundException {
		AbstractSequence<AminoAcidCompound> seq = super.getSequence(sequence.toUpperCase(Locale.ENGLISH), index);
		seq.setUserCollection(getStringCase(sequence));
		return seq;
	}


	/**
	 * Assumes all compounds were uppercase
	 * @see org.biojava.nbio.core.sequence.io.ProteinSequenceCreator#getSequence(java.util.List)
	 */
	@Override
	public AbstractSequence<AminoAcidCompound> getSequence(
			List<AminoAcidCompound> list) {
		AbstractSequence<AminoAcidCompound> seq =super.getSequence(list);
		Collection<Object> strCase = new ArrayList<Object>(seq.getLength());
		for(int i=0;i<seq.getLength();i++) {
			strCase.add(true);
		}
		seq.setUserCollection(strCase);
		return seq;
	}

	/**
	 * Returns a list of Booleans of the same length as the input, specifying
	 * whether each character was uppercase or not.
	 * @param str A string. Should not contain unicode supplemental characters.
	 * @return a list of Booleans of the same length as the input, specifying
	 * whether each character was uppercase or not.
	 * This list contains only Booleans.
	 */
	private static List<Object> getStringCase(String str) {
		List<Object> types = new ArrayList<Object>(str.length());
		for(int i=0;i<str.length();i++) {
			types.add(Character.isUpperCase(str.charAt(i)));
		}
		return types;
	}

	public static void main(String[] args) throws Exception {
		CasePreservingProteinSequenceCreator creator = new CasePreservingProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet());
		AbstractSequence<AminoAcidCompound> seq = creator.getSequence("aaAA",0);
		logger.info("Sequence: {}", seq.getSequenceAsString()); //"AAAA"
		logger.info("User Collection: {}", seq.getUserCollection()); //"[false, false, true, true]"
	}

	/**
	 * Takes a {@link ProteinSequence} which was created by a
	 * {@link CasePreservingProteinSequenceCreator}. Uses the case info
	 * stored in the user collection to modify the output array.
	 *
	 * <p>Sets elements of the output array which correspond to lowercase letters
	 * to null.
	 *
	 * @param seq Input sequence with case stored as the user collection
	 * @param out
	 */
	public static void setLowercaseToNull( ProteinSequence seq,
			Object[] out) {
		// should have been set by seq creator
		Collection<Object> userCollection = seq.getUserCollection();
		if(userCollection == null)
			throw new IllegalArgumentException("Sequence doesn't contain valid case info");
		if(userCollection.size() != out.length)
			throw new IllegalArgumentException("Sequence length doesn't math output array length");

		int pos = 0;
		for(Object isAligned : userCollection) {
			assert(isAligned instanceof Boolean);
			if(!(Boolean)isAligned) {
				out[pos] = null;
			}
			pos++;
		}
	}
}
