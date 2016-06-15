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
 *
 */
package org.biojava.nbio.structure.io.mmcif.model;


/**
 *  Data items in the ENTITY_SRC_GEN category record details of
 *  the source from which the entity was obtained in cases
 *  where the source was genetically manipulated.  The
 *  following are treated separately:  items pertaining to the tissue
 *  from which the gene was obtained, items pertaining to the host
 *  organism for gene expression and items pertaining to the actual
 *  producing organism (plasmid).
 *
 * @author Andreas Prlic
 *
 */
public class EntitySrcGen {
	String  entity_id;
	String  expression_system_id;
	String  gene_src_common_name;
	String  gene_src_details	;
	String  gene_src_dev_stage	 ;
	String  gene_src_genus	 ;
	String  gene_src_species	;
	String  gene_src_strain	 ;
	String  gene_src_tissue	 ;
	String  gene_src_tissue_fraction;
	String  host_org_common_name	 ;
	String  host_org_details	 ;
	String  host_org_genus	 ;
	String  host_org_species;
	String  host_org_strain	 ;
	String  pdbx_src_id;
	String  pdbx_seq_type;
	String  pdbx_alt_source_flag;
	String  pdbx_beg_seq_num;
	String  pdbx_end_seq_num;
	String  pdbx_description;
	String  pdbx_gene_src_atcc;
	String  pdbx_gene_src_cell	;
	String  pdbx_gene_src_cell_line;
	String  pdbx_gene_src_cellular_location;
	String  pdbx_gene_src_fragment	 ;
	String  pdbx_gene_src_gene	 ;
	String  pdbx_gene_src_ncbi_taxonomy_id;
	String  pdbx_gene_src_organ	 ;
	String  pdbx_gene_src_organelle ;
	String  pdbx_gene_src_plasmid	 ;
	String  pdbx_gene_src_plasmid_name	 ;
	String  pdbx_gene_src_scientific_name;
	String  pdbx_gene_src_variant	 ;
	String  pdbx_host_org_atcc	 ;
	String  pdbx_host_org_cell	 ;
	String  pdbx_host_org_cell_line ;
	String  pdbx_host_org_cellular_location ;
	String  pdbx_host_org_culture_collection ;
	String  pdbx_host_org_gene	  ;
	String  pdbx_host_org_ncbi_taxonomy_id ;
	String  pdbx_host_org_organ	 ;
	String  pdbx_host_org_organelle ;
	String  pdbx_host_org_scientific_name ;
	String  pdbx_host_org_strain	  ;
	String  pdbx_host_org_tissue	 ;
	String  pdbx_host_org_tissue_fraction ;
	String  pdbx_host_org_variant	 ;
	String  pdbx_host_org_vector	 ;
	String  pdbx_host_org_vector_type;
	String  plasmid_details	 ;
	String  plasmid_name	 ;
	String  start_construct_id ;
	public String getEntity_id() {
		return entity_id;
	}
	public void setEntity_id(String entity_id) {
		this.entity_id = entity_id;
	}
	public String getExpression_system_id() {
		return expression_system_id;
	}
	public void setExpression_system_id(String expression_system_id) {
		this.expression_system_id = expression_system_id;
	}
	public String getGene_src_common_name() {
		return gene_src_common_name;
	}
	public void setGene_src_common_name(String gene_src_common_name) {
		this.gene_src_common_name = gene_src_common_name;
	}
	public String getGene_src_details() {
		return gene_src_details;
	}
	public void setGene_src_details(String gene_src_details) {
		this.gene_src_details = gene_src_details;
	}
	public String getGene_src_dev_stage() {
		return gene_src_dev_stage;
	}
	public void setGene_src_dev_stage(String gene_src_dev_stage) {
		this.gene_src_dev_stage = gene_src_dev_stage;
	}
	public String getGene_src_genus() {
		return gene_src_genus;
	}
	public void setGene_src_genus(String gene_src_genus) {
		this.gene_src_genus = gene_src_genus;
	}
	public String getGene_src_species() {
		return gene_src_species;
	}
	public void setGene_src_species(String gene_src_species) {
		this.gene_src_species = gene_src_species;
	}
	public String getGene_src_strain() {
		return gene_src_strain;
	}
	public void setGene_src_strain(String gene_src_strain) {
		this.gene_src_strain = gene_src_strain;
	}
	public String getGene_src_tissue() {
		return gene_src_tissue;
	}
	public void setGene_src_tissue(String gene_src_tissue) {
		this.gene_src_tissue = gene_src_tissue;
	}
	public String getGene_src_tissue_fraction() {
		return gene_src_tissue_fraction;
	}
	public void setGene_src_tissue_fraction(String gene_src_tissue_fraction) {
		this.gene_src_tissue_fraction = gene_src_tissue_fraction;
	}
	public String getHost_org_common_name() {
		return host_org_common_name;
	}
	public void setHost_org_common_name(String host_org_common_name) {
		this.host_org_common_name = host_org_common_name;
	}
	public String getHost_org_details() {
		return host_org_details;
	}
	public void setHost_org_details(String host_org_details) {
		this.host_org_details = host_org_details;
	}
	public String getHost_org_genus() {
		return host_org_genus;
	}
	public void setHost_org_genus(String host_org_genus) {
		this.host_org_genus = host_org_genus;
	}
	public String getHost_org_species() {
		return host_org_species;
	}
	public void setHost_org_species(String host_org_species) {
		this.host_org_species = host_org_species;
	}
	public String getHost_org_strain() {
		return host_org_strain;
	}
	public void setHost_org_strain(String host_org_strain) {
		this.host_org_strain = host_org_strain;
	}
	public String getPdbx_src_id() {
		return pdbx_src_id;
	}
	public void setPdbx_src_id(String pdbx_src_id) {
		this.pdbx_src_id = pdbx_src_id;
	}
	public String getPdbx_seq_type() {
		return pdbx_seq_type;
	}
	public void setPdbx_seq_type(String pdbx_seq_type) {
		this.pdbx_seq_type = pdbx_seq_type;
	}
	/**
	 * @return the pdbx_alt_source_flag
	 */
	public String getPdbx_alt_source_flag() {
		return pdbx_alt_source_flag;
	}
	/**
	 * @param pdbx_alt_source_flag the pdbx_alt_source_flag to set
	 */
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
	public String getPdbx_description() {
		return pdbx_description;
	}
	public void setPdbx_description(String pdbx_description) {
		this.pdbx_description = pdbx_description;
	}
	public String getPdbx_gene_src_atcc() {
		return pdbx_gene_src_atcc;
	}
	public void setPdbx_gene_src_atcc(String pdbx_gene_src_atcc) {
		this.pdbx_gene_src_atcc = pdbx_gene_src_atcc;
	}
	public String getPdbx_gene_src_cell() {
		return pdbx_gene_src_cell;
	}
	public void setPdbx_gene_src_cell(String pdbx_gene_src_cell) {
		this.pdbx_gene_src_cell = pdbx_gene_src_cell;
	}
	public String getPdbx_gene_src_cell_line() {
		return pdbx_gene_src_cell_line;
	}
	public void setPdbx_gene_src_cell_line(String pdbx_gene_src_cell_line) {
		this.pdbx_gene_src_cell_line = pdbx_gene_src_cell_line;
	}
	public String getPdbx_gene_src_cellular_location() {
		return pdbx_gene_src_cellular_location;
	}
	public void setPdbx_gene_src_cellular_location(
			String pdbx_gene_src_cellular_location) {
		this.pdbx_gene_src_cellular_location = pdbx_gene_src_cellular_location;
	}
	public String getPdbx_gene_src_fragment() {
		return pdbx_gene_src_fragment;
	}
	public void setPdbx_gene_src_fragment(String pdbx_gene_src_fragment) {
		this.pdbx_gene_src_fragment = pdbx_gene_src_fragment;
	}
	public String getPdbx_gene_src_gene() {
		return pdbx_gene_src_gene;
	}
	public void setPdbx_gene_src_gene(String pdbx_gene_src_gene) {
		this.pdbx_gene_src_gene = pdbx_gene_src_gene;
	}
	public String getPdbx_gene_src_ncbi_taxonomy_id() {
		return pdbx_gene_src_ncbi_taxonomy_id;
	}
	public void setPdbx_gene_src_ncbi_taxonomy_id(
			String pdbx_gene_src_ncbi_taxonomy_id) {
		this.pdbx_gene_src_ncbi_taxonomy_id = pdbx_gene_src_ncbi_taxonomy_id;
	}
	public String getPdbx_gene_src_organ() {
		return pdbx_gene_src_organ;
	}
	public void setPdbx_gene_src_organ(String pdbx_gene_src_organ) {
		this.pdbx_gene_src_organ = pdbx_gene_src_organ;
	}
	public String getPdbx_gene_src_organelle() {
		return pdbx_gene_src_organelle;
	}
	public void setPdbx_gene_src_organelle(String pdbx_gene_src_organelle) {
		this.pdbx_gene_src_organelle = pdbx_gene_src_organelle;
	}
	public String getPdbx_gene_src_plasmid() {
		return pdbx_gene_src_plasmid;
	}
	public void setPdbx_gene_src_plasmid(String pdbx_gene_src_plasmid) {
		this.pdbx_gene_src_plasmid = pdbx_gene_src_plasmid;
	}
	public String getPdbx_gene_src_plasmid_name() {
		return pdbx_gene_src_plasmid_name;
	}
	public void setPdbx_gene_src_plasmid_name(String pdbx_gene_src_plasmid_name) {
		this.pdbx_gene_src_plasmid_name = pdbx_gene_src_plasmid_name;
	}
	public String getPdbx_gene_src_scientific_name() {
		return pdbx_gene_src_scientific_name;
	}
	public void setPdbx_gene_src_scientific_name(
			String pdbx_gene_src_scientific_name) {
		this.pdbx_gene_src_scientific_name = pdbx_gene_src_scientific_name;
	}
	public String getPdbx_gene_src_variant() {
		return pdbx_gene_src_variant;
	}
	public void setPdbx_gene_src_variant(String pdbx_gene_src_variant) {
		this.pdbx_gene_src_variant = pdbx_gene_src_variant;
	}
	public String getPdbx_host_org_atcc() {
		return pdbx_host_org_atcc;
	}
	public void setPdbx_host_org_atcc(String pdbx_host_org_atcc) {
		this.pdbx_host_org_atcc = pdbx_host_org_atcc;
	}
	public String getPdbx_host_org_cell() {
		return pdbx_host_org_cell;
	}
	public void setPdbx_host_org_cell(String pdbx_host_org_cell) {
		this.pdbx_host_org_cell = pdbx_host_org_cell;
	}
	public String getPdbx_host_org_cell_line() {
		return pdbx_host_org_cell_line;
	}
	public void setPdbx_host_org_cell_line(String pdbx_host_org_cell_line) {
		this.pdbx_host_org_cell_line = pdbx_host_org_cell_line;
	}
	public String getPdbx_host_org_cellular_location() {
		return pdbx_host_org_cellular_location;
	}
	public void setPdbx_host_org_cellular_location(
			String pdbx_host_org_cellular_location) {
		this.pdbx_host_org_cellular_location = pdbx_host_org_cellular_location;
	}
	public String getPdbx_host_org_culture_collection() {
		return pdbx_host_org_culture_collection;
	}
	public void setPdbx_host_org_culture_collection(
			String pdbx_host_org_culture_collection) {
		this.pdbx_host_org_culture_collection = pdbx_host_org_culture_collection;
	}
	public String getPdbx_host_org_gene() {
		return pdbx_host_org_gene;
	}
	public void setPdbx_host_org_gene(String pdbx_host_org_gene) {
		this.pdbx_host_org_gene = pdbx_host_org_gene;
	}
	public String getPdbx_host_org_ncbi_taxonomy_id() {
		return pdbx_host_org_ncbi_taxonomy_id;
	}
	public void setPdbx_host_org_ncbi_taxonomy_id(
			String pdbx_host_org_ncbi_taxonomy_id) {
		this.pdbx_host_org_ncbi_taxonomy_id = pdbx_host_org_ncbi_taxonomy_id;
	}
	public String getPdbx_host_org_organ() {
		return pdbx_host_org_organ;
	}
	public void setPdbx_host_org_organ(String pdbx_host_org_organ) {
		this.pdbx_host_org_organ = pdbx_host_org_organ;
	}
	public String getPdbx_host_org_organelle() {
		return pdbx_host_org_organelle;
	}
	public void setPdbx_host_org_organelle(String pdbx_host_org_organelle) {
		this.pdbx_host_org_organelle = pdbx_host_org_organelle;
	}
	public String getPdbx_host_org_scientific_name() {
		return pdbx_host_org_scientific_name;
	}
	public void setPdbx_host_org_scientific_name(
			String pdbx_host_org_scientific_name) {
		this.pdbx_host_org_scientific_name = pdbx_host_org_scientific_name;
	}
	public String getPdbx_host_org_strain() {
		return pdbx_host_org_strain;
	}
	public void setPdbx_host_org_strain(String pdbx_host_org_strain) {
		this.pdbx_host_org_strain = pdbx_host_org_strain;
	}
	public String getPdbx_host_org_tissue() {
		return pdbx_host_org_tissue;
	}
	public void setPdbx_host_org_tissue(String pdbx_host_org_tissue) {
		this.pdbx_host_org_tissue = pdbx_host_org_tissue;
	}
	public String getPdbx_host_org_tissue_fraction() {
		return pdbx_host_org_tissue_fraction;
	}
	public void setPdbx_host_org_tissue_fraction(
			String pdbx_host_org_tissue_fraction) {
		this.pdbx_host_org_tissue_fraction = pdbx_host_org_tissue_fraction;
	}
	public String getPdbx_host_org_variant() {
		return pdbx_host_org_variant;
	}
	public void setPdbx_host_org_variant(String pdbx_host_org_variant) {
		this.pdbx_host_org_variant = pdbx_host_org_variant;
	}
	public String getPdbx_host_org_vector() {
		return pdbx_host_org_vector;
	}
	public void setPdbx_host_org_vector(String pdbx_host_org_vector) {
		this.pdbx_host_org_vector = pdbx_host_org_vector;
	}
	public String getPdbx_host_org_vector_type() {
		return pdbx_host_org_vector_type;
	}
	public void setPdbx_host_org_vector_type(String pdbx_host_org_vector_type) {
		this.pdbx_host_org_vector_type = pdbx_host_org_vector_type;
	}
	public String getPlasmid_details() {
		return plasmid_details;
	}
	public void setPlasmid_details(String plasmid_details) {
		this.plasmid_details = plasmid_details;
	}
	public String getPlasmid_name() {
		return plasmid_name;
	}
	public void setPlasmid_name(String plasmid_name) {
		this.plasmid_name = plasmid_name;
	}
	public String getStart_construct_id() {
		return start_construct_id;
	}
	public void setStart_construct_id(String start_construct_id) {
		this.start_construct_id = start_construct_id;
	}




}
