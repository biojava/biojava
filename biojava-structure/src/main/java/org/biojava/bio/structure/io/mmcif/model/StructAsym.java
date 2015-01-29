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
 * created at Jun 1, 2008
 */
package org.biojava.bio.structure.io.mmcif.model;

/** Contains the data for _struct_asym 
 * 
 * @author Andreas Prlic
 * @since 1.7
 *
 */
public class StructAsym extends AbstractBean{
	String id; 
	String pdbx_blank_PDB_chainid_flag; 
	String pdbx_modified; 
	String entity_id; 
	String details;
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public String getPdbx_blank_PDB_chainid_flag() {
		return pdbx_blank_PDB_chainid_flag;
	}
	public void setPdbx_blank_PDB_chainid_flag(String pdbx_blank_PDB_chainid_flag) {
		this.pdbx_blank_PDB_chainid_flag = pdbx_blank_PDB_chainid_flag;
	}
	public String getPdbx_modified() {
		return pdbx_modified;
	}
	public void setPdbx_modified(String pdbx_modified) {
		this.pdbx_modified = pdbx_modified;
	}
	public String getEntity_id() {
		return entity_id;
	}
	public void setEntity_id(String entity_id) {
		this.entity_id = entity_id;
	}
	public String getDetails() {
		return details;
	}
	public void setDetails(String details) {
		this.details = details;
	}
	
	
}
