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


public class Symmetry extends AbstractBean {

	String entry_id;
	@CIFLabel(label="space_group_name_H-M")
	String space_group_name_H_M;
	@CIFLabel(label="pdbx_full_space_group_name_H-M")
	String pdbx_full_space_group_name_H_M;
	String cell_setting;
	String Int_Tables_number;
	String space_group_name_Hall;
	public String getEntry_id() {
		return entry_id;
	}
	public void setEntry_id(String entry_id) {
		this.entry_id = entry_id;
	}
	public String getSpace_group_name_H_M() {
		return space_group_name_H_M;
	}
	public void setSpace_group_name_H_M(String space_group_name_H_M) {
		this.space_group_name_H_M = space_group_name_H_M;
	}
	public String getPdbx_full_space_group_name_H_M() {
		return pdbx_full_space_group_name_H_M;
	}
	public void setPdbx_full_space_group_name_H_M(
			String pdbx_full_space_group_name_H_M) {
		this.pdbx_full_space_group_name_H_M = pdbx_full_space_group_name_H_M;
	}
	public String getCell_setting() {
		return cell_setting;
	}
	public void setCell_setting(String cell_setting) {
		this.cell_setting = cell_setting;
	}
	public String getInt_Tables_number() {
		return Int_Tables_number;
	}
	public void setInt_Tables_number(String int_Tables_number) {
		Int_Tables_number = int_Tables_number;
	}
	public String getSpace_group_name_Hall() {
		return space_group_name_Hall;
	}
	public void setSpace_group_name_Hall(String space_group_name_Hall) {
		this.space_group_name_Hall = space_group_name_Hall;
	}
}
