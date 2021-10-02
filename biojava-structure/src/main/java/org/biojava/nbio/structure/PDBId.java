package org.biojava.nbio.structure;

import java.io.Serializable;
import java.util.regex.Pattern;

public class PDBId implements Comparable<PDBId>, Serializable{
	
	private static final long serialVersionUID = -5400143283145477754L;
	public enum Behavior{
		PREFER_SHORT,
		PREFER_EXTENDED
	}
	private static Behavior defaultShorteningBehaviour = Behavior.PREFER_SHORT;

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

	public static final Pattern PATTERN_SHORT_PDBID = Pattern.compile("[1-9]\\p{Alnum}{3}");
//	public static final Pattern PATTERN_EXTENDED_PDBID = Pattern.compile("(PDB|pdb)_\\d{4}[1-9]\\p{Alnum}{3}");//Shall we allow lower case?
	public static final Pattern PATTERN_EXTENDED_PDBID = Pattern.compile("PDB_\\d{4}[1-9]\\p{Alnum}{3}"); 

	/**
	 * Keeps the ID in UPPER CASE, in a reduced form (without the <code>PDB_</code> prefix).
	 */
	private String idCode;
	private static final String XXXX_STRING = "XXXX";

	public static final PDBId XXXX = getXxxx();

	/*
	 * using a static method instead of static initializer because the later cant't
	 * easily handle initializing a final field with an object whose constructor
	 * throws a checked exception correctly.
	 */	
	private static PDBId getXxxx(){ 
		try {
			return new PDBId(XXXX_STRING);
		} catch (PDBIdException e) {
			return null; //will never happen
		}
	}
	
	public PDBId(String id) throws PDBIdException {
		if (id == null) {
			throw new PDBIdException("ID can not be null");
		}
		if(accept_xxxx && XXXX_STRING.equalsIgnoreCase(id)) {// the only exception
			this.idCode = toInternalFormat(XXXX_STRING);
		}else {
			this.idCode = toInternalFormat(id);
		}
	}
	
	public static boolean isShortPDBID(String id) throws NullPointerException {
		return PATTERN_SHORT_PDBID.matcher(id).matches();
	}
	
	public static boolean isExtendedPDBID(String id) throws NullPointerException {
		return PATTERN_EXTENDED_PDBID.matcher(id).matches();
	}
	
	/**Checks whether an Extended PDBId is shortable, <i>assuming it is a valid extended PDBId</i>.
	 * If you are not sure the String represents a valid extended PDBId, use {@link #isExtendedPDBID(String)} first.
	 * @see #isExtendedPDBID(String)
	 * @param extendedId
	 * @return
	 */
	public static boolean isShortCompatible(String extendedId) {
		return extendedId.length() >= 8 && extendedId.substring(0, 8).equals/*IgnoreCase*/("PDB_0000");
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
		return this.getId().equals(((PDBId)obj).getId());
	}
	
	@Override
	protected Object clone() throws CloneNotSupportedException {
		try {
			return new PDBId(this.getId());
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
	 * @param b when it equals <code>Behavior.PREFER_SHORT</code>, the class will try to produce the short ID whenever possible.
	 * @return The PDBId in short format if possible and <code>b</code> equals <code>Behavior.PREFER_SHORT</code>, the extended PdBID form otherwise.
	 */
	public String getId(Behavior b) {
		if (b == Behavior.PREFER_SHORT && isInternalShortCompatible(idCode))
			return internalToShortNoCheck(idCode);
		return internalToExtendedFormat(idCode);
	}
	
	public String getShortId() throws PDBIdException{
		if(isInternalShortCompatible(idCode)) {
			return internalToShortNoCheck(idCode);
		} else {
			throw new PDBIdException("ID ("+getId()+") is not short format compatible");
		}
	}
	
	public static String toExtendedId(String shortId) throws PDBIdException{
		if (isShortPDBID(shortId) || XXXX_STRING.equalsIgnoreCase(shortId)) {
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

	private static boolean isInternalShortCompatible(String intId) {
		return intId.substring(0, 4).equalsIgnoreCase("0000");
	}
	
	private static String toInternalFormat(String id) throws PDBIdException {
		if (isShortPDBID(id) || XXXX_STRING.equalsIgnoreCase(id)) {
			return ("0000"+id).toUpperCase();
		}else if (isExtendedPDBID(id)) {
			return id.substring(4).toUpperCase();
		} else {
			throw new PDBIdException("Unknown format ["+id+"]");
		}
	}
	
	private static String internalToExtendedFormat(String id){
		return "PDB_" + id;
	}
	
	private static String internalToShortNoCheck(String extendedId) {
		return extendedId.substring(4).toUpperCase();
	}
	
	private static String toShortNoCheck(String extendedId) {
		return extendedId.substring(8).toUpperCase();
	}

	@Override
	public int compareTo(PDBId o) {
		return this.idCode.compareTo(o.idCode);
	}
	
}
