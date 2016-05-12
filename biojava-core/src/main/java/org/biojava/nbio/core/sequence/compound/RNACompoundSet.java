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
package org.biojava.nbio.core.sequence.compound;

import org.biojava.nbio.core.sequence.template.AbstractNucleotideCompoundSet;

/**
 *
 * @author Andy Yates
 */
public class RNACompoundSet extends AbstractNucleotideCompoundSet<NucleotideCompound> {

	private static class InitaliseOnDemand {
		public static final RNACompoundSet INSTANCE = new RNACompoundSet();
	}
	public static RNACompoundSet getRNACompoundSet() {
		return InitaliseOnDemand.INSTANCE;
	}

	public RNACompoundSet() {
		addNucleotideCompound("A", "U");
		addNucleotideCompound("U", "A");
		addNucleotideCompound("G", "C");
		addNucleotideCompound("C", "G");
		addNucleotideCompound("N", "N");
		addNucleotideCompound("-", "-");
	}

	@Override
public NucleotideCompound newNucleotideCompound(String base, String complement, String... equivalents) {
		return new NucleotideCompound(base, this, complement);
	}
}
