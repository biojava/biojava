/**
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
 * Created on Dec 7, 2013
 * Created by Douglas Myers-Turnbull
 *
 */
package org.biojava.nbio.structure.io.sifts;

/**
 * An entry in the chain-level SIFTS mapping between UniProt and the PDB.
 * @author dmyersturnbull
 * @see SiftsChainToUniprotMapping
 * @since 3.0.7
 * @see SiftsResidue Which is a distinct concept
 */
public class SiftsChainEntry {

	private final String chainId;
	private final String pdbEnd;
	private final String pdbId;
	private final String pdbStart;
	private final String seqresEnd;
	private final String seqresStart;
	private final String uniprotEnd;
	private final String uniProtId;
	private final String uniprotStart;

	public SiftsChainEntry(String pdbId, String chainId, String uniProtId, String seqresStart, String seqresEnd,
			String pdbStart, String pdbEnd, String uniprotStart, String uniprotEnd) {
		super();
		this.pdbId = pdbId;
		this.chainId = chainId;
		this.uniProtId = uniProtId;
		this.seqresStart = seqresStart;
		this.seqresEnd = seqresEnd;
		this.pdbStart = pdbStart;
		this.pdbEnd = pdbEnd;
		this.uniprotStart = uniprotStart;
		this.uniprotEnd = uniprotEnd;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		SiftsChainEntry other = (SiftsChainEntry) obj;
		if (chainId == null) {
			if (other.chainId != null) return false;
		} else if (!chainId.equals(other.chainId)) return false;
		if (pdbEnd == null) {
			if (other.pdbEnd != null) return false;
		} else if (!pdbEnd.equals(other.pdbEnd)) return false;
		if (pdbId == null) {
			if (other.pdbId != null) return false;
		} else if (!pdbId.equals(other.pdbId)) return false;
		if (pdbStart == null) {
			if (other.pdbStart != null) return false;
		} else if (!pdbStart.equals(other.pdbStart)) return false;
		if (seqresEnd == null) {
			if (other.seqresEnd != null) return false;
		} else if (!seqresEnd.equals(other.seqresEnd)) return false;
		if (seqresStart == null) {
			if (other.seqresStart != null) return false;
		} else if (!seqresStart.equals(other.seqresStart)) return false;
		if (uniProtId == null) {
			if (other.uniProtId != null) return false;
		} else if (!uniProtId.equals(other.uniProtId)) return false;
		if (uniprotEnd == null) {
			if (other.uniprotEnd != null) return false;
		} else if (!uniprotEnd.equals(other.uniprotEnd)) return false;
		if (uniprotStart == null) {
			if (other.uniprotStart != null) return false;
		} else if (!uniprotStart.equals(other.uniprotStart)) return false;
		return true;
	}

	public String getChainId() {
		return chainId;
	}

	/**
	 * @return A residue number
	 */
	public String getPdbEnd() {
		return pdbEnd;
	}

	public String getPdbId() {
		return pdbId;
	}

	/**
	 * @return A residue number
	 */
	public String getPdbStart() {
		return pdbStart;
	}

	/**
	 * @return A residue number
	 */
	public String getSeqresEnd() {
		return seqresEnd;
	}

	/**
	 * @return A residue number
	 */
	public String getSeqresStart() {
		return seqresStart;
	}

	public String getUniprotEnd() {
		return uniprotEnd;
	}

	public String getUniProtId() {
		return uniProtId;
	}

	public String getUniprotStart() {
		return uniprotStart;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + (chainId == null ? 0 : chainId.hashCode());
		result = prime * result + (pdbEnd == null ? 0 : pdbEnd.hashCode());
		result = prime * result + (pdbId == null ? 0 : pdbId.hashCode());
		result = prime * result + (pdbStart == null ? 0 : pdbStart.hashCode());
		result = prime * result + (seqresEnd == null ? 0 : seqresEnd.hashCode());
		result = prime * result + (seqresStart == null ? 0 : seqresStart.hashCode());
		result = prime * result + (uniProtId == null ? 0 : uniProtId.hashCode());
		result = prime * result + (uniprotEnd == null ? 0 : uniprotEnd.hashCode());
		result = prime * result + (uniprotStart == null ? 0 : uniprotStart.hashCode());
		return result;
	}

	@Override
	public String toString() {
		return "SiftsChainToUniprotEntry [pdbId=" + pdbId + ", chainName=" + chainId + ", uniProtId=" + uniProtId
				+ ", seqresStart=" + seqresStart + ", seqresEnd=" + seqresEnd + ", pdbStart=" + pdbStart + ", pdbEnd="
				+ pdbEnd + ", uniprotStart=" + uniprotStart + ", uniprotEnd=" + uniprotEnd + "]";
	}

}
