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
package org.biojava.nbio.structure.contact;

import java.io.Serializable;

public class AtomIdentifier implements Serializable {

	private static final long serialVersionUID = 1L;

	private int pdbSerial;
	private String chainId;

	public AtomIdentifier(int pdbSerial, String chainId) {
		this.pdbSerial = pdbSerial;
		this.chainId = chainId;
	}

	public int getPdbSerial() {
		return pdbSerial;
	}

	public void setPdbSerial(int pdbSerial) {
		this.pdbSerial = pdbSerial;
	}

	public String getChainId() {
		return chainId;
	}

	public void setChainId(String chainId) {
		this.chainId = chainId;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((chainId == null) ? 0 : chainId.hashCode());
		result = prime * result + pdbSerial;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		AtomIdentifier other = (AtomIdentifier) obj;
		if (chainId == null) {
			if (other.chainId != null)
				return false;
		} else if (!chainId.equals(other.chainId))
			return false;
		if (pdbSerial != other.pdbSerial)
			return false;
		return true;
	}

	@Override
	public String toString() {
		return " [" + pdbSerial + " - "
				+ chainId + "]";
	}


}
