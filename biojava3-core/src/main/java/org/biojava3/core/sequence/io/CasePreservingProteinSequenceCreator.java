package org.biojava3.core.sequence.io;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Locale;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;

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

	public CasePreservingProteinSequenceCreator(
			CompoundSet<AminoAcidCompound> compoundSet) {
		super(compoundSet);
	}

	/**
	 * 
	 * @see org.biojava3.core.sequence.io.ProteinSequenceCreator#getSequence(org.biojava3.core.sequence.template.ProxySequenceReader, long)
	 */
	@Override
	public AbstractSequence<AminoAcidCompound> getSequence(
			ProxySequenceReader<AminoAcidCompound> proxyLoader, long index) {
		AbstractSequence<AminoAcidCompound> seq = super.getSequence(proxyLoader, index);
		seq.setUserCollection(getStringCase(proxyLoader.getSequenceAsString()));
		return seq;
	}

	/* (non-Javadoc)
	 * @see org.biojava3.core.sequence.io.ProteinSequenceCreator#getSequence(java.lang.String, long)
	 */
	@Override
	public AbstractSequence<AminoAcidCompound> getSequence(String sequence,
			long index) {
		AbstractSequence<AminoAcidCompound> seq = super.getSequence(sequence.toUpperCase(Locale.ENGLISH), index);
		seq.setUserCollection(getStringCase(sequence));
		return seq;
	}
	

	/**
	 * Assumes all compounds were uppercase
	 * @see org.biojava3.core.sequence.io.ProteinSequenceCreator#getSequence(java.util.List)
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
	
	public static void main(String[] args) {
		CasePreservingProteinSequenceCreator creator = new CasePreservingProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet());
		AbstractSequence<AminoAcidCompound> seq = creator.getSequence("aaAA",0);
		System.out.println(seq.getSequenceAsString()); //"AAAA"
		System.out.println(seq.getUserCollection()); //"[false, false, true, true]"
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
