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
package org.biojava.nbio.structure.io.mmcif.model;

/**
 * A class containing the _atom_sites data. Equivalent to the SCALE records in PDB files.
 * 
 * 
 * @author Jose Duarte
 *
 */
public class AtomSites extends AbstractBean {
	
	String entry_id;
	
	// to my knowledge this field is not used - JD 2016-11-22
	String Cartn_transform_axes;
	
	@CIFLabel(label="fract_transf_matrix[1][1]")
	String fract_transf_matrix11;

	@CIFLabel(label="fract_transf_matrix[1][2]")
	String fract_transf_matrix12;

	@CIFLabel(label="fract_transf_matrix[1][3]")
	String fract_transf_matrix13;


	@CIFLabel(label="fract_transf_matrix[2][1]")
	String fract_transf_matrix21;

	@CIFLabel(label="fract_transf_matrix[2][2]")
	String fract_transf_matrix22;

	@CIFLabel(label="fract_transf_matrix[2][3]")
	String fract_transf_matrix23;
	

	@CIFLabel(label="fract_transf_matrix[3][1]")
	String fract_transf_matrix31;

	@CIFLabel(label="fract_transf_matrix[3][2]")
	String fract_transf_matrix32;

	@CIFLabel(label="fract_transf_matrix[3][3]")
	String fract_transf_matrix33;


	@CIFLabel(label="fract_transf_vector[1]")
	String fract_transf_vector1;
	
	@CIFLabel(label="fract_transf_vector[2]")
	String fract_transf_vector2;

	@CIFLabel(label="fract_transf_vector[3]")
	String fract_transf_vector3;
	
	// these fields are unusual but appear in some entries like 5e5j - JD 2016-11-22
	@CIFLabel(label="Cartn_transf_matrix[1][1]")
	String Cartn_transf_matrix11;

	@CIFLabel(label="Cartn_transf_matrix[1][2]")
	String Cartn_transf_matrix12;

	@CIFLabel(label="Cartn_transf_matrix[1][3]")
	String Cartn_transf_matrix13;

	@CIFLabel(label="Cartn_transf_matrix[2][1]")
	String Cartn_transf_matrix21;

	@CIFLabel(label="Cartn_transf_matrix[2][2]")
	String Cartn_transf_matrix22;

	@CIFLabel(label="Cartn_transf_matrix[2][3]")
	String Cartn_transf_matrix23;

	@CIFLabel(label="Cartn_transf_matrix[3][1]")
	String Cartn_transf_matrix31;

	@CIFLabel(label="Cartn_transf_matrix[3][2]")
	String Cartn_transf_matrix32;

	@CIFLabel(label="Cartn_transf_matrix[3][3]")
	String Cartn_transf_matrix33;
	
	@CIFLabel(label="Cartn_transf_vector[1]")
	String Cartn_transf_vector1;

	@CIFLabel(label="Cartn_transf_vector[2]")
	String Cartn_transf_vector2;

	@CIFLabel(label="Cartn_transf_vector[3]")
	String Cartn_transf_vector3;

	
	public String getEntry_id() {
		return entry_id;
	}
	public void setEntry_id(String entry_id) {
		this.entry_id = entry_id;
	}
	/**
	 * @return the cartn_transform_axes
	 */
	public String getCartn_transform_axes() {
		return Cartn_transform_axes;
	}
	/**
	 * @param cartn_transform_axes the cartn_transform_axes to set
	 */
	public void setCartn_transform_axes(String cartn_transform_axes) {
		Cartn_transform_axes = cartn_transform_axes;
	}
	/**
	 * @return the fract_transf_matrix11
	 */
	public String getFract_transf_matrix11() {
		return fract_transf_matrix11;
	}
	/**
	 * @param fract_transf_matrix11 the fract_transf_matrix11 to set
	 */
	public void setFract_transf_matrix11(String fract_transf_matrix11) {
		this.fract_transf_matrix11 = fract_transf_matrix11;
	}
	/**
	 * @return the fract_transf_matrix12
	 */
	public String getFract_transf_matrix12() {
		return fract_transf_matrix12;
	}
	/**
	 * @param fract_transf_matrix12 the fract_transf_matrix12 to set
	 */
	public void setFract_transf_matrix12(String fract_transf_matrix12) {
		this.fract_transf_matrix12 = fract_transf_matrix12;
	}
	/**
	 * @return the fract_transf_matrix13
	 */
	public String getFract_transf_matrix13() {
		return fract_transf_matrix13;
	}
	/**
	 * @param fract_transf_matrix13 the fract_transf_matrix13 to set
	 */
	public void setFract_transf_matrix13(String fract_transf_matrix13) {
		this.fract_transf_matrix13 = fract_transf_matrix13;
	}
	/**
	 * @return the fract_transf_matrix21
	 */
	public String getFract_transf_matrix21() {
		return fract_transf_matrix21;
	}
	/**
	 * @param fract_transf_matrix21 the fract_transf_matrix21 to set
	 */
	public void setFract_transf_matrix21(String fract_transf_matrix21) {
		this.fract_transf_matrix21 = fract_transf_matrix21;
	}
	/**
	 * @return the fract_transf_matrix22
	 */
	public String getFract_transf_matrix22() {
		return fract_transf_matrix22;
	}
	/**
	 * @param fract_transf_matrix22 the fract_transf_matrix22 to set
	 */
	public void setFract_transf_matrix22(String fract_transf_matrix22) {
		this.fract_transf_matrix22 = fract_transf_matrix22;
	}
	/**
	 * @return the fract_transf_matrix23
	 */
	public String getFract_transf_matrix23() {
		return fract_transf_matrix23;
	}
	/**
	 * @param fract_transf_matrix23 the fract_transf_matrix23 to set
	 */
	public void setFract_transf_matrix23(String fract_transf_matrix23) {
		this.fract_transf_matrix23 = fract_transf_matrix23;
	}
	/**
	 * @return the fract_transf_matrix31
	 */
	public String getFract_transf_matrix31() {
		return fract_transf_matrix31;
	}
	/**
	 * @param fract_transf_matrix31 the fract_transf_matrix31 to set
	 */
	public void setFract_transf_matrix31(String fract_transf_matrix31) {
		this.fract_transf_matrix31 = fract_transf_matrix31;
	}
	/**
	 * @return the fract_transf_matrix32
	 */
	public String getFract_transf_matrix32() {
		return fract_transf_matrix32;
	}
	/**
	 * @param fract_transf_matrix32 the fract_transf_matrix32 to set
	 */
	public void setFract_transf_matrix32(String fract_transf_matrix32) {
		this.fract_transf_matrix32 = fract_transf_matrix32;
	}
	/**
	 * @return the fract_transf_matrix33
	 */
	public String getFract_transf_matrix33() {
		return fract_transf_matrix33;
	}
	/**
	 * @param fract_transf_matrix33 the fract_transf_matrix33 to set
	 */
	public void setFract_transf_matrix33(String fract_transf_matrix33) {
		this.fract_transf_matrix33 = fract_transf_matrix33;
	}
	/**
	 * @return the fract_transf_vector1
	 */
	public String getFract_transf_vector1() {
		return fract_transf_vector1;
	}
	/**
	 * @param fract_transf_vector1 the fract_transf_vector1 to set
	 */
	public void setFract_transf_vector1(String fract_transf_vector1) {
		this.fract_transf_vector1 = fract_transf_vector1;
	}
	/**
	 * @return the fract_transf_vector2
	 */
	public String getFract_transf_vector2() {
		return fract_transf_vector2;
	}
	/**
	 * @param fract_transf_vector2 the fract_transf_vector2 to set
	 */
	public void setFract_transf_vector2(String fract_transf_vector2) {
		this.fract_transf_vector2 = fract_transf_vector2;
	}
	/**
	 * @return the fract_transf_vector3
	 */
	public String getFract_transf_vector3() {
		return fract_transf_vector3;
	}
	/**
	 * @param fract_transf_vector3 the fract_transf_vector3 to set
	 */
	public void setFract_transf_vector3(String fract_transf_vector3) {
		this.fract_transf_vector3 = fract_transf_vector3;
	}
	/**
	 * @return the cartn_transf_matrix11
	 */
	public String getCartn_transf_matrix11() {
		return Cartn_transf_matrix11;
	}
	/**
	 * @param cartn_transf_matrix11 the cartn_transf_matrix11 to set
	 */
	public void setCartn_transf_matrix11(String cartn_transf_matrix11) {
		Cartn_transf_matrix11 = cartn_transf_matrix11;
	}
	/**
	 * @return the cartn_transf_matrix12
	 */
	public String getCartn_transf_matrix12() {
		return Cartn_transf_matrix12;
	}
	/**
	 * @param cartn_transf_matrix12 the cartn_transf_matrix12 to set
	 */
	public void setCartn_transf_matrix12(String cartn_transf_matrix12) {
		Cartn_transf_matrix12 = cartn_transf_matrix12;
	}
	/**
	 * @return the cartn_transf_matrix13
	 */
	public String getCartn_transf_matrix13() {
		return Cartn_transf_matrix13;
	}
	/**
	 * @param cartn_transf_matrix13 the cartn_transf_matrix13 to set
	 */
	public void setCartn_transf_matrix13(String cartn_transf_matrix13) {
		Cartn_transf_matrix13 = cartn_transf_matrix13;
	}
	/**
	 * @return the cartn_transf_matrix21
	 */
	public String getCartn_transf_matrix21() {
		return Cartn_transf_matrix21;
	}
	/**
	 * @param cartn_transf_matrix21 the cartn_transf_matrix21 to set
	 */
	public void setCartn_transf_matrix21(String cartn_transf_matrix21) {
		Cartn_transf_matrix21 = cartn_transf_matrix21;
	}
	/**
	 * @return the cartn_transf_matrix22
	 */
	public String getCartn_transf_matrix22() {
		return Cartn_transf_matrix22;
	}
	/**
	 * @param cartn_transf_matrix22 the cartn_transf_matrix22 to set
	 */
	public void setCartn_transf_matrix22(String cartn_transf_matrix22) {
		Cartn_transf_matrix22 = cartn_transf_matrix22;
	}
	/**
	 * @return the cartn_transf_matrix23
	 */
	public String getCartn_transf_matrix23() {
		return Cartn_transf_matrix23;
	}
	/**
	 * @param cartn_transf_matrix23 the cartn_transf_matrix23 to set
	 */
	public void setCartn_transf_matrix23(String cartn_transf_matrix23) {
		Cartn_transf_matrix23 = cartn_transf_matrix23;
	}
	/**
	 * @return the cartn_transf_matrix31
	 */
	public String getCartn_transf_matrix31() {
		return Cartn_transf_matrix31;
	}
	/**
	 * @param cartn_transf_matrix31 the cartn_transf_matrix31 to set
	 */
	public void setCartn_transf_matrix31(String cartn_transf_matrix31) {
		Cartn_transf_matrix31 = cartn_transf_matrix31;
	}
	/**
	 * @return the cartn_transf_matrix32
	 */
	public String getCartn_transf_matrix32() {
		return Cartn_transf_matrix32;
	}
	/**
	 * @param cartn_transf_matrix32 the cartn_transf_matrix32 to set
	 */
	public void setCartn_transf_matrix32(String cartn_transf_matrix32) {
		Cartn_transf_matrix32 = cartn_transf_matrix32;
	}
	/**
	 * @return the cartn_transf_matrix33
	 */
	public String getCartn_transf_matrix33() {
		return Cartn_transf_matrix33;
	}
	/**
	 * @param cartn_transf_matrix33 the cartn_transf_matrix33 to set
	 */
	public void setCartn_transf_matrix33(String cartn_transf_matrix33) {
		Cartn_transf_matrix33 = cartn_transf_matrix33;
	}
	/**
	 * @return the cartn_transf_vector1
	 */
	public String getCartn_transf_vector1() {
		return Cartn_transf_vector1;
	}
	/**
	 * @param cartn_transf_vector1 the cartn_transf_vector1 to set
	 */
	public void setCartn_transf_vector1(String cartn_transf_vector1) {
		Cartn_transf_vector1 = cartn_transf_vector1;
	}
	/**
	 * @return the cartn_transf_vector2
	 */
	public String getCartn_transf_vector2() {
		return Cartn_transf_vector2;
	}
	/**
	 * @param cartn_transf_vector2 the cartn_transf_vector2 to set
	 */
	public void setCartn_transf_vector2(String cartn_transf_vector2) {
		Cartn_transf_vector2 = cartn_transf_vector2;
	}
	/**
	 * @return the cartn_transf_vector3
	 */
	public String getCartn_transf_vector3() {
		return Cartn_transf_vector3;
	}
	/**
	 * @param cartn_transf_vector3 the cartn_transf_vector3 to set
	 */
	public void setCartn_transf_vector3(String cartn_transf_vector3) {
		Cartn_transf_vector3 = cartn_transf_vector3;
	}

}
