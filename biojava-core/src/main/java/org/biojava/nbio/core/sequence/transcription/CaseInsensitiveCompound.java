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

import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.Compound;

/**
 * Attempts to wrap compounds so it is possible to view them
 * in a case insensitive manner
 */
public class CaseInsensitiveCompound implements Compound {

	private final NucleotideCompound compound;

	public CaseInsensitiveCompound(NucleotideCompound compound) {
		this.compound = compound;
	}

	@Override
public boolean equalsIgnoreCase(Compound compound) {
		if (compound == null) {
			return false;
		}
		if (!(compound instanceof CaseInsensitiveCompound)) {
			return false;
		}
		CaseInsensitiveCompound them = (CaseInsensitiveCompound) compound;
		return toString().equalsIgnoreCase(them.toString());
	}

	@Override
public boolean equals(Object obj) {
		if (obj == null) {
			return false;
		}
		if (!(obj instanceof CaseInsensitiveCompound)) {
			return false;
		}
		return equalsIgnoreCase((Compound)obj);
	}

	@Override
public int hashCode() {
		return toString().toUpperCase().hashCode();
	}

	public NucleotideCompound getUnderlyingCompound() {
		return this.compound;
	}

	@Override
public String getDescription() {
		return getUnderlyingCompound().getDescription();
	}

	@Override
public String getLongName() {
		return getUnderlyingCompound().getLongName();
	}

	@Override
public Float getMolecularWeight() {
		return getUnderlyingCompound().getMolecularWeight();
	}

	@Override
public String getShortName() {
		return getUnderlyingCompound().getShortName();
	}

	@Override
public String toString() {
		return getUnderlyingCompound().toString();
	}

	@Override
public void setDescription(String description) {
		//Nothing
	}

	@Override
public void setLongName(String longName) {
		//Nothing
	}

	@Override
public void setMolecularWeight(Float molecularWeight) {
		//Nothing
	}

	@Override
public void setShortName(String shortName) {
		//Nothing
	}
}
