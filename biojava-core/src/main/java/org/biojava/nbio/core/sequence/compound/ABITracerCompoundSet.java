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
 * Created on 05-06-2018
 */

package org.biojava.nbio.core.sequence.compound;

import org.biojava.nbio.core.sequence.template.AbstractNucleotideCompoundSet;

/**
 * @author Maximilian Greil
 */
public class ABITracerCompoundSet extends AbstractNucleotideCompoundSet<NucleotideCompound> {

	private static class InitaliseOnDemand {
		public static final ABITracerCompoundSet INSTANCE = new ABITracerCompoundSet();
	}

	public static ABITracerCompoundSet getABITracerCompoundSet() {
		return InitaliseOnDemand.INSTANCE;
	}

	public ABITracerCompoundSet() {
		addNucleotideCompound("A", "T");
		addNucleotideCompound("T", "A");
		addNucleotideCompound("G", "C");
		addNucleotideCompound("C", "G");
		addNucleotideCompound("N", "N");
		addNucleotideCompound("K", "K");
		addNucleotideCompound("Y", "Y");
		addNucleotideCompound("R", "R");
		addNucleotideCompound("-", "-");
	}

	@Override
	public NucleotideCompound newNucleotideCompound(String base, String complement, String... equivalents) {
		if(equivalents.length == 0) {
			return new NucleotideCompound(base, this, complement);
		}
		else {
			NucleotideCompound[] compounds = new NucleotideCompound[equivalents.length];
			for(int i=0; i<compounds.length; i++) {
				compounds[i] = getCompoundForString(equivalents[i]);
			}
			return new NucleotideCompound(base, this, complement, compounds);
		}
	}
}
