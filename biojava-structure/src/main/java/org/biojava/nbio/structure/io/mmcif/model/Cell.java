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

public class Cell extends AbstractBean {

	String entry_id;
	String length_a;
	String length_b;
	String length_c;
	String angle_alpha;
	String angle_beta;
	String angle_gamma;
	String Z_PDB;
	String pdbx_unique_axis;

	// some PDB entries like 1aac have the extra esd fields
	String length_a_esd;
	String length_b_esd;
	String length_c_esd;
	String angle_alpha_esd;
	String angle_beta_esd;
	String angle_gamma_esd;
	
	String volume;

	public String getEntry_id() {
		return entry_id;
	}
	public void setEntry_id(String entry_id) {
		this.entry_id = entry_id;
	}
	public String getLength_a() {
		return length_a;
	}
	public void setLength_a(String length_a) {
		this.length_a = length_a;
	}
	public String getLength_b() {
		return length_b;
	}
	public void setLength_b(String length_b) {
		this.length_b = length_b;
	}
	public String getLength_c() {
		return length_c;
	}
	public void setLength_c(String length_c) {
		this.length_c = length_c;
	}
	public String getAngle_alpha() {
		return angle_alpha;
	}
	public void setAngle_alpha(String angle_alpha) {
		this.angle_alpha = angle_alpha;
	}
	public String getAngle_beta() {
		return angle_beta;
	}
	public void setAngle_beta(String angle_beta) {
		this.angle_beta = angle_beta;
	}
	public String getAngle_gamma() {
		return angle_gamma;
	}
	public void setAngle_gamma(String angle_gamma) {
		this.angle_gamma = angle_gamma;
	}
	public String getZ_PDB() {
		return Z_PDB;
	}
	public void setZ_PDB(String z_PDB) {
		Z_PDB = z_PDB;
	}
	public String getPdbx_unique_axis() {
		return pdbx_unique_axis;
	}
	public void setPdbx_unique_axis(String pdbx_unique_axis) {
		this.pdbx_unique_axis = pdbx_unique_axis;
	}
	public String getLength_a_esd() {
		return length_a_esd;
	}
	public void setLength_a_esd(String length_a_esd) {
		this.length_a_esd = length_a_esd;
	}
	public String getLength_b_esd() {
		return length_b_esd;
	}
	public void setLength_b_esd(String length_b_esd) {
		this.length_b_esd = length_b_esd;
	}
	public String getLength_c_esd() {
		return length_c_esd;
	}
	public void setLength_c_esd(String length_c_esd) {
		this.length_c_esd = length_c_esd;
	}
	public String getAngle_alpha_esd() {
		return angle_alpha_esd;
	}
	public void setAngle_alpha_esd(String angle_alpha_esd) {
		this.angle_alpha_esd = angle_alpha_esd;
	}
	public String getAngle_beta_esd() {
		return angle_beta_esd;
	}
	public void setAngle_beta_esd(String angle_beta_esd) {
		this.angle_beta_esd = angle_beta_esd;
	}
	public String getAngle_gamma_esd() {
		return angle_gamma_esd;
	}
	public void setAngle_gamma_esd(String angle_gamma_esd) {
		this.angle_gamma_esd = angle_gamma_esd;
	}
	public String getVolume() {
		return volume;
	}
	public void setVolume(String volume) {
		this.volume = volume;
	}


}
