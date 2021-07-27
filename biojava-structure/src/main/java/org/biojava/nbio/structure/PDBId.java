package org.biojava.nbio.structure;

import java.io.Serializable;

public class PDBId implements Comparable<PDBId>, Serializable{
	
	private static final String XXXX = "XXXX";
	private static final long serialVersionUID = -5400143283145477754L;
	public static final boolean PREFER_SHORT = true;
	public static final boolean PREFER_EXTENDED = false;
	private static boolean defaultShorteningBehaviour = PREFER_SHORT;

	private static boolean accept_xxxx = true;

	public static class PDBIdException extends StructureException{
		private static final long serialVersionUID = -9166852283492713918L;
		public PDBIdException(String message) {
			super(message);
		}
		public PDBIdException(String message, Throwable cause) {
			super(message, cause);
		}
		public PDBIdException(Throwable cause) {
			super(cause);
		}
	}

	public static final String SHORT_PDBID_FORMRMAT		= "\\d{1}\\p{Alnum}{3}";
	public static final String EXTENDED_PDBID_FORMRMAT	= "PDB_\\d{5}\\p{Alnum}{3}";

	/**
	 * keeps the ID in UPPER CASE.
	 */
	private String idCode;
	
	public PDBId(String id) throws PDBIdException {
		if (id == null) {
			throw new PDBIdException("ID can not be null");
		}
		if(accept_xxxx && XXXX.equalsIgnoreCase(id)) {// the only exception
			this.idCode = toExtendedId(XXXX);
		}else if (isExtendedPDBID(id)) {
			this.idCode = id.toUpperCase();
		} else {
			this.idCode = toExtendedId(id);
		}
	}
	
	public static boolean isShortPDBID(String id) throws NullPointerException {
		return id.matches(SHORT_PDBID_FORMRMAT);
	}
	
	public static boolean isExtendedPDBID(String id) throws NullPointerException {
		return id.matches(EXTENDED_PDBID_FORMRMAT);
	}
	
	public static boolean isShortCompatible(String extendedId) {
		return extendedId.substring(0, 8).equalsIgnoreCase("PDB_0000");
	}
	
	@Override
	public int hashCode() {
		return idCode.hashCode();
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		// We are sure they are both objects of the same class and their ID is in the same (UPPER) case.
		return this.hashCode() == obj.hashCode();
	}
	
	@Override
	protected Object clone() throws CloneNotSupportedException {
		try {
			return new PDBId(idCode);
		} catch (PDBIdException e) {
			e.printStackTrace();  // Can never happen. Do nothing
			throw new CloneNotSupportedException("Unexpected error");
		}
	}

	@Override
	public String toString() {
		return getId();
	}

	/**By default this function will try to get the PDBId in the short (4 letters) format.
	 * If not possible, it will return the long format.
	 * N.B. This default behavior may change later;
	 * @return the PDBId code, preferably in short format.
	 */
	public String getId() {
		return getId(defaultShorteningBehaviour);
	}
	
	/**
	 * @param prefereShort When <code>true</code>, the class will try to produce the short ID whenever possible.
	 * @return The PDBId in short format if possible and <code>prefereShort</code> is <code>true</code>, the extended PdBID form otherwise.
	 */
	public String getId(boolean prefereShort) {
		if (prefereShort && isShortCompatible(idCode))
			return toShortNoCheck(idCode);
		return idCode;
	}
	
	public String getShortId() throws PDBIdException{
		return toShortId(idCode);
	}
	
	public static String toExtendedId(String shortId) throws PDBIdException{
		if (isShortPDBID(shortId) || XXXX.equalsIgnoreCase(shortId)) {
			return ("PDB_0000"+shortId).toUpperCase();
		}else if (isExtendedPDBID(shortId)) {
			return shortId.toUpperCase();
		} else {
			throw new PDBIdException("Unknown format ["+shortId+"]");
		}
	}
	
	public static String toShortId(String extendedId) throws PDBIdException{
		if (isExtendedPDBID(extendedId) && isShortCompatible(extendedId)) {
			return toShortNoCheck(extendedId);
		} else if (isShortPDBID(extendedId)) {
			return extendedId.toUpperCase();
		} else {
			throw new PDBIdException("Conversion not possible of ID ["+extendedId+"]");
		}
	}

	private static String toShortNoCheck(String extendedId) {
		return extendedId.substring(8).toUpperCase();
	}
	

	@Override
	public int compareTo(PDBId o) {
		return this.idCode.compareTo(o.idCode);
	}
	
}
