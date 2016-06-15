/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do  t have a copy,
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
 * created at Aug 12, 2013
 * Author: Andreas Prlic
 */

package org.biojava.nbio.structure.io.mmcif.model;

/** Data items in the ENTITY_SRC_NAT category record details of
			   the source from which the entity was obtained in cases
			   where the entity was isolated directly from a natural tissue.
 */
public class EntitySrcNat {
	String  common_name	;
	String  details	;
	String  entity_id ;
	String  genus	 ;
	String  pdbx_atcc	;
	String  pdbx_cell	 ;
	String  pdbx_cell_line ;
	String  pdbx_cellular_location;
	String  pdbx_fragment	 ;
	String  pdbx_ncbi_taxonomy_id;
	String  pdbx_organ	 ;
	String  pdbx_organelle;
	String  pdbx_organism_scientific;
	String  pdbx_plasmid_details	 ;
	String  pdbx_plasmid_name	 ;
	String  pdbx_secretion	 ;
	String  pdbx_variant	 ;
	String  pdbx_src_id;
	String  pdbx_alt_source_flag;
	String  pdbx_beg_seq_num;
	String  pdbx_end_seq_num;
	String  pdbx_leaving_atom_flag;
	String  species	 ;
	String  strain	 ;
	String  tissue	 ;
	String  tissue_fraction;

	public String getCommon_name() {
		return common_name;
	}
	public void setCommon_name(String common_name) {
		this.common_name = common_name;
	}
	public String getDetails() {
		return details;
	}
	public void setDetails(String details) {
		this.details = details;
	}
	public String getEntity_id() {
		return entity_id;
	}
	public void setEntity_id(String entity_id) {
		this.entity_id = entity_id;
	}
	public String getGenus() {
		return genus;
	}
	public void setGenus(String genus) {
		this.genus = genus;
	}
	public String getPdbx_atcc() {
		return pdbx_atcc;
	}
	public void setPdbx_atcc(String pdbx_atcc) {
		this.pdbx_atcc = pdbx_atcc;
	}
	public String getPdbx_cell() {
		return pdbx_cell;
	}
	public void setPdbx_cell(String pdbx_cell) {
		this.pdbx_cell = pdbx_cell;
	}
	public String getPdbx_cell_line() {
		return pdbx_cell_line;
	}
	public void setPdbx_cell_line(String pdbx_cell_line) {
		this.pdbx_cell_line = pdbx_cell_line;
	}
	public String getPdbx_cellular_location() {
		return pdbx_cellular_location;
	}
	public void setPdbx_cellular_location(String pdbx_cellular_location) {
		this.pdbx_cellular_location = pdbx_cellular_location;
	}
	public String getPdbx_fragment() {
		return pdbx_fragment;
	}
	public void setPdbx_fragment(String pdbx_fragment) {
		this.pdbx_fragment = pdbx_fragment;
	}
	public String getPdbx_ncbi_taxonomy_id() {
		return pdbx_ncbi_taxonomy_id;
	}
	public void setPdbx_ncbi_taxonomy_id(String pdbx_ncbi_taxonomy_id) {
		this.pdbx_ncbi_taxonomy_id = pdbx_ncbi_taxonomy_id;
	}
	public String getPdbx_organ() {
		return pdbx_organ;
	}
	public void setPdbx_organ(String pdbx_organ) {
		this.pdbx_organ = pdbx_organ;
	}
	public String getPdbx_organelle() {
		return pdbx_organelle;
	}
	public void setPdbx_organelle(String pdbx_organelle) {
		this.pdbx_organelle = pdbx_organelle;
	}
	public String getPdbx_organism_scientific() {
		return pdbx_organism_scientific;
	}
	public void setPdbx_organism_scientific(String pdbx_organism_scientific) {
		this.pdbx_organism_scientific = pdbx_organism_scientific;
	}
	public String getPdbx_plasmid_details() {
		return pdbx_plasmid_details;
	}
	public void setPdbx_plasmid_details(String pdbx_plasmid_details) {
		this.pdbx_plasmid_details = pdbx_plasmid_details;
	}
	public String getPdbx_plasmid_name() {
		return pdbx_plasmid_name;
	}
	public void setPdbx_plasmid_name(String pdbx_plasmid_name) {
		this.pdbx_plasmid_name = pdbx_plasmid_name;
	}
	public String getPdbx_secretion() {
		return pdbx_secretion;
	}
	public void setPdbx_secretion(String pdbx_secretion) {
		this.pdbx_secretion = pdbx_secretion;
	}
	public String getPdbx_variant() {
		return pdbx_variant;
	}
	public void setPdbx_variant(String pdbx_variant) {
		this.pdbx_variant = pdbx_variant;
	}
	public String getSpecies() {
		return species;
	}
	public void setSpecies(String species) {
		this.species = species;
	}
	public String getStrain() {
		return strain;
	}
	public void setStrain(String strain) {
		this.strain = strain;
	}
	public String getTissue() {
		return tissue;
	}
	public void setTissue(String tissue) {
		this.tissue = tissue;
	}
	public String getTissue_fraction() {
		return tissue_fraction;
	}
	public void setTissue_fraction(String tissue_fraction) {
		this.tissue_fraction = tissue_fraction;
	}

	public String getPdbx_src_id() {
		return pdbx_src_id;
	}

	public void setPdbx_src_id(String pdbx_src_id) {
		this.pdbx_src_id = pdbx_src_id;
	}

	public String getPdbx_alt_source_flag() {
		return pdbx_alt_source_flag;
	}

	public void setPdbx_alt_source_flag(String pdbx_alt_source_flag) {
		this.pdbx_alt_source_flag = pdbx_alt_source_flag;
	}

	public String getPdbx_beg_seq_num() {
		return pdbx_beg_seq_num;
	}

	public void setPdbx_beg_seq_num(String pdbx_beg_seq_num) {
		this.pdbx_beg_seq_num = pdbx_beg_seq_num;
	}

	public String getPdbx_end_seq_num() {
		return pdbx_end_seq_num;
	}

	public void setPdbx_end_seq_num(String pdbx_end_seq_num) {
		this.pdbx_end_seq_num = pdbx_end_seq_num;
	}

	public String getPdbx_leaving_atom_flag() {
		return pdbx_leaving_atom_flag;
	}

	public void setPdbx_leaving_atom_flag(String pdbx_leaving_atom_flag) {
		this.pdbx_leaving_atom_flag = pdbx_leaving_atom_flag;
	}
}
