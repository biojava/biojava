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
package org.biojava.nbio.structure.symmetry.misc;

import java.util.List;

public class ChainSignature implements Comparable<ChainSignature> {
	private int chainCount = 0;
	private String representative = "";
	private String compositionId = "";
	private List<String> chainIds = null;

	public ChainSignature(String representative, int chainCount, List<String> chainIds) {
		this.representative = representative;
		this.chainCount = chainCount;
		this.chainIds = chainIds;
	}

	public String getRepresentative() {
		return representative;
	}

	public List<String> getChainIds() {
		return chainIds;
	}

	public String getCompositionId() {
		return compositionId;
	}

	public void setCompositionId(String compositionId) {
		this.compositionId = compositionId;
	}

	@Override
	public boolean equals (Object obj) {
		if (this == obj) {
			return true;
		}
		if((obj == null) || (obj.getClass() != this.getClass())) return false;

		ChainSignature other = (ChainSignature) obj;
		if (representative == null) {
			return false;
		}
		return chainCount == other.chainCount && representative.equals(other.representative);
	}

	@Override
	public int compareTo(ChainSignature other) {
		if (other.chainCount < this.chainCount) {
			return -1;
		}
		if (other.chainCount > this.chainCount) {
			return 1;
		}
		return this.representative.compareTo(other.representative);
	}

	@Override
	public int hashCode() {
		int hash = 7;
		hash = 31 * hash + chainCount;
		hash = 31 * hash + (representative == null ? 0 : representative.hashCode());
		return hash;
	}

	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		builder.append("(");
		builder.append(representative);
		builder.append(")");
		if (chainCount> 1) {
			builder.append(chainCount);
		}
		return builder.toString();
	}

}
