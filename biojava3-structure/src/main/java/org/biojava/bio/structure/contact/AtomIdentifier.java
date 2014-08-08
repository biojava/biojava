package org.biojava.bio.structure.contact;

public class AtomIdentifier {

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
