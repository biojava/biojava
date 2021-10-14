package org.biojava.nbio.structure;

import java.io.Serializable;
import java.util.regex.Pattern;

public class PdbId implements Comparable<PdbId>, Serializable{
	
	private static final String PREFIX_PDB_ = "PDB_";
	private static final String STRING_0000 = "0000";
	private static final String PDB_0000 = PREFIX_PDB_ + STRING_0000;
	private static final String XXXX_STRING = "XXXX";

	private static final long serialVersionUID = -5400143283145477754L;
	public enum Behavior{
		PREFER_SHORT,
		PREFER_EXTENDED
	}
	private static Behavior defaultShorteningBehaviour = Behavior.PREFER_SHORT;

	private static boolean accept_xxxx = true;


	public static final Pattern PATTERN_SHORT_PDBID = Pattern.compile("[1-9]\\p{Alnum}{3}");
//	public static final Pattern PATTERN_EXTENDED_PDBID = Pattern.compile("(PDB|pdb)_\\d{4}[1-9]\\p{Alnum}{3}");//Shall we allow lower case?
	public static final Pattern PATTERN_EXTENDED_PDBID = Pattern.compile("PDB_\\d{4}[1-9]\\p{Alnum}{3}"); 

	/**
	 * Keeps the ID in UPPER CASE, in a reduced form (without the <code>PDB_</code> prefix).
	 */
	private String idCode;

	public static final PdbId XXXX = new PdbId(XXXX_STRING);

	public PdbId(String id){
		if (id == null) {
			throw new IllegalArgumentException("ID can not be null");
		}
		if(accept_xxxx && XXXX_STRING.equalsIgnoreCase(id)) {// the only exception
			this.idCode = toInternalFormat(XXXX_STRING);
		}else {
			this.idCode = toInternalFormat(id);
		}
	}
	
	/**
	 * @param id
	 * @return
	 * @throws NullPointerException
	 */
	public static boolean isValidShortPdbId(String id) {
		return PATTERN_SHORT_PDBID.matcher(id).matches();
	}
	
	/**
	 * @param id
	 * @return
	 * @throws NullPointerException
	 */
	public static boolean isValidExtendedPdbId(String id) {
		return PATTERN_EXTENDED_PDBID.matcher(id).matches();
	}
	
	/**Checks whether an Extended PDB ID is shortable, <i>assuming it is a valid extended PDB ID</i>.
	 * If you are not sure the String represents a valid extended PdbId, use {@link #isValidExtendedPdbId(String)} first.
	 * @see #isValidExtendedPdbId(String)
	 * @param extendedId
	 * @return
	 */
	public static boolean isShortCompatible(String extendedId) {
		return extendedId.length() >= 8 && extendedId.substring(0, 8).equals/*IgnoreCase*/(PDB_0000);
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
		return this.getId().equals(((PdbId)obj).getId());
	}
	
	@Override
	protected Object clone() throws CloneNotSupportedException {
		return new PdbId(this.getId());
	}

	@Override
	public String toString() {
		return getId();
	}

	/**By default this function will try to get the PdbId in the short (4 letters) format.
	 * If not possible, it will return the long format.
	 * N.B. This default behavior may change later;
	 * @return the PdbId code, preferably in short format.
	 */
	public String getId() {
		return getId(defaultShorteningBehaviour);
	}
	
	/**
	 * @param b when it equals <code>Behavior.PREFER_SHORT</code>, the class will try to produce the short ID whenever possible.
	 * @return The PdbId in short format if possible and <code>b</code> equals <code>Behavior.PREFER_SHORT</code>, the extended PdBID form otherwise.
	 */
	public String getId(Behavior b) {
		if (b == Behavior.PREFER_SHORT && isInternalShortCompatible(idCode))
			return internalToShortNoCheck(idCode);
		return PREFIX_PDB_ + idCode;
	}
	
	public String getShortId() throws StructureException{
		if(isInternalShortCompatible(idCode)) {
			return internalToShortNoCheck(idCode);
		} else {
			throw new StructureException("ID (" + getId() + ") is not short format compatible");
		}
	}
	
	public static String toExtendedId(String shortId) throws StructureException{
		if (isValidShortPdbId(shortId) || XXXX_STRING.equalsIgnoreCase(shortId)) {
			return PDB_0000 + shortId.toUpperCase();
		}else if (isValidExtendedPdbId(shortId)) {
			return shortId.toUpperCase();
		} else {
			throw new StructureException("Unknown format ["+shortId+"]");
		}
	}
	
	public static String toShortId(String extendedId) throws StructureException{
		if (isValidExtendedPdbId(extendedId) && isShortCompatible(extendedId)) {
			return extendedId.substring(8).toUpperCase();
		} else if (isValidShortPdbId(extendedId)) {
			return extendedId.toUpperCase();
		} else {
			throw new StructureException("Conversion not possible of ID ["+extendedId+"]");
		}
	}

	private static boolean isInternalShortCompatible(String intId) {
		return intId.substring(0, 4).equals(STRING_0000);
	}
	
	private static String toInternalFormat(String id) throws IllegalArgumentException {
		if (isValidShortPdbId(id) || XXXX_STRING.equalsIgnoreCase(id)) {
			return STRING_0000  + id.toUpperCase();
		}else if (isValidExtendedPdbId(id)) {
			return id.substring(4).toUpperCase();
		} else {
			throw new IllegalArgumentException("Unknown format [" + id + "]");
		}
	}
	
	private static String internalToShortNoCheck(String extendedId) {
		return extendedId.substring(4).toUpperCase();
	}
	
	@Override
	public int compareTo(PdbId o) {
		return this.idCode.compareTo(o.idCode);
	}
	
}
