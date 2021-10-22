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
package org.biojava.nbio.structure;

import java.io.Serializable;
import java.util.regex.Pattern;

/**
 * A wrapper class for the PDB identifier.
 * 
 * It handles conversion between current (short) <code>[1-9][0-9A-Z]{3}</code> and
 * upcoming (extended) <code>PDB_\d{4}[1-9][09-A-Z]</code> PDB ID format.<br>
 * Instances of this class are <em>immutable</em>.<br>
 * Creation of PdBId instance follows strict PDB ID convention.
 * There is only one exception to this rule which is <b>XXXX</b>. XXXX objects 
 * are not considered equal (unless they are the one and the same object).
 * @author Amr ALHOSSARY
 * @since 6.0.0
 *
 */
public class PdbId implements Comparable<PdbId>, Serializable{
	
	private static final String PREFIX_PDB_ = "PDB_";
	private static final String STRING_0000 = "0000";
	private static final String PDB_0000 = PREFIX_PDB_ + STRING_0000;
	public static final String XXXX_STRING = "XXXX";
	private static final long serialVersionUID = -5577433541029161356L;

	/**
	 * How the PDB ID output/conversion should go, if possible.
	 */
	public enum Behavior{
		/**
		 * Try to produce short PDB ID. If failed, produce extended PDB ID.
		 */
		PREFER_SHORT,
		/**
		 * Always produce Extended PDB ID.
		 */
		PREFER_EXTENDED
	}
	private static Behavior defaultShorteningBehaviour = Behavior.PREFER_SHORT;

	private static boolean accept_xxxx = true;


	/**
	 * A regular expression that matches a PDB ID in the short format.
	 */
	public static final Pattern PATTERN_SHORT_PDBID = Pattern.compile("[1-9]\\p{Alnum}{3}");
	/**
	 * A regular expression that matches a PDB ID in the extended format.
	 */
	public static final Pattern PATTERN_EXTENDED_PDBID = Pattern.compile("PDB_\\d{4}[1-9]\\p{Alnum}{3}"); 

	/**
	 * Keeps the ID in <b>UPPER CASE</b>, in a <em>reduced</em> form (without the <code>PDB_</code> prefix).
	 */
	private String idCode;

	private static final String XXXX_INTERNAL = toInternalFormat(XXXX_STRING);
	/**
	 * @param id A <i>valid</i> PDB ID in either <i>short (case insensitive)</i> or <i>extended</i> format.
	 * @throws IllegalArgumentException If <code>id</code> is not a valid identifier.
	 * @throws NullPointerException  If <code>id</code> is <code>null</code>.
	 */
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
	 * Check whether <code>id</code> represents a valid PDB ID in the <em>short</em> format.
	 * @param id Prospect ID
	 * @return <code>true</code> if <code>id</code> is a valid short PDB ID, <code>false</code> otherwise.
	 * @throws NullPointerException if <code>id</code> is <code>null</code>.
	 * @see #isValidExtendedPdbId(String)
	 */
	public static boolean isValidShortPdbId(String id) {
		return PATTERN_SHORT_PDBID.matcher(id).matches();
	}
	
	/**
	 * Check whether <code>id</code> represents a valid PDB ID in the <em>extended</em> format.
	 * @param id Prospect ID
	 * @return <code>true</code> if <code>id</code> is a valid extended PDB ID, <code>false</code> otherwise.
	 * @throws NullPointerException if <code>id</code> is <code>null</code>.
	 * @see #isValidShortPdbId(String)
	 */
	public static boolean isValidExtendedPdbId(String id) {
		return PATTERN_EXTENDED_PDBID.matcher(id).matches();
	}
	
	/**
	 * Checks whether an Extended PDB ID is shortable, <i>assuming it is a valid extended PDB ID</i>.
	 * If you are not sure the String represents a valid extended PdbId, use {@link #isValidExtendedPdbId(String)} first.
	 * @see #isValidExtendedPdbId(String)
	 * @param extendedId the supposedly valid extended PDB ID.
	 * @return <code>true</code> if <code>extendedId</code> starts with "PDB_0000", <code>false</code> otherwise.
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
		// Time for checking for XXXX. If this is XXXX or obj is XXXX return false.
		// Checking whether this is XXXX only is enough because Checking whether
		//   obj is XXXX is implicitly included in the next check. So no need to
		//   waste some machine cycles checking it here.
		if(XXXX_INTERNAL.equals(this.idCode))
			return false;
		// We are sure they are both objects of the same class and their respective IDs are in the same (UPPER) case.
		return this.idCode.equals(((PdbId)obj).idCode);
	}
	
	@Override
	protected Object clone() throws CloneNotSupportedException {
		return new PdbId(this.getId());
	}

	@Override
	public String toString() {
		return getId();
	}

	/**
	 * Get a <code>String</code> representation of this PdbId instance.<br>
	 * By default this function will try to get the PdbId in the short (4 letters) format.
	 * If not possible, it will return the long format.
	 * N.B. This default behavior may change later;
	 * @return the PdbId code, preferably in short format.
	 */
	public String getId() {
		return getId(defaultShorteningBehaviour);
	}
	
	/**
	 * Get a <code>String</code> representation of this PdbId instance, using the <i>passed in</i> behavior.<br>
	 * @param b when it equals <code>Behavior.PREFER_SHORT</code>, the class will try to produce the short ID whenever possible.
	 * @return The PdbId in short format if possible and <code>b</code> equals <code>Behavior.PREFER_SHORT</code>, the extended PDB ID form otherwise.
	 */
	public String getId(Behavior b) {
		if (b == Behavior.PREFER_SHORT && isInternalShortCompatible(idCode))
			return internalToShortNoCheck(idCode);
		return PREFIX_PDB_ + idCode;
	}
	
	/**
	 * Get the PDB Id in the sort format. Throws an exception if the conversion is not possible.<br>
	 * Use this method only if you know that this PDB ID is shortable.
	 * @return the PDB ID in the short format.
	 * @throws StructureException if the conversion was not possible.
	 */
	public String getShortId() throws StructureException{
		if(isInternalShortCompatible(idCode)) {
			return internalToShortNoCheck(idCode);
		} else {
			throw new StructureException("ID (" + getId() + ") is not short format compatible");
		}
	}
	
	/**
	 * Converts <code>shortId</code> to the PDB ID extended format.
	 * If <code>shortId</code> is a valid short PDB ID (or XXXX), it would be converted to an extended ID,
	 * if <code>shortId</code> is a valid extended PDB ID, it would be returned in UPPER CASE,
	 * a {@link StructureException} is thrown otherwise.
	 * @param shortId the PDB ID to convert to extended format
	 * @return the ID in the extended UPPER CASE format.
	 * @throws StructureException if the conversion was not possible.
	 */
	public static String toExtendedId(String shortId) throws StructureException{
		if (isValidShortPdbId(shortId) || XXXX_STRING.equalsIgnoreCase(shortId)) {
			return PDB_0000 + shortId.toUpperCase();
		}else if (isValidExtendedPdbId(shortId)) {
			return shortId.toUpperCase();
		} else {
			throw new StructureException("Unknown format ["+shortId+"]");
		}
	}
	
	/**
	 * Converts <code>extendedId</code> to the PDB ID short format.
	 * If <code>extendedId</code> is a valid extended PDB ID, it would be converted to a short ID,
	 * if <code>extendedId</code> is a valid short PDB ID, it would be returned in UPPER CASE,
	 * a {@link StructureException} is thrown otherwise.
	 * @param extendedId the PDB ID to convert to short format
	 * @return the ID in the short UPPER CASE format.
	 * @throws StructureException if the conversion was not possible.
	 */
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
		//We know that both idCode fields are 8 UPPER CASE characters strings.
		return this.idCode.compareTo(o.idCode);
	}
	
}
