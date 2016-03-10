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
 * created at Aug 12, 2013
 * Author: Andreas Prlic
 */

package org.biojava.nbio.structure.io.mmcif.model;

/**
 * PDBX_ENTITY_SRC_SYN records the details about each chemically
 * synthesized molecule (entity) in the asymmetric unit.
 * @author Andreas Prlic
 *
 */
public class EntitySrcSyn {
	String  details;
	String  entity_id;
	String  ncbi_taxonomy_id;
	String  organism_common_name;
	String  organism_scientific;
	String  strain;
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
	public String getNcbi_taxonomy_id() {
		return ncbi_taxonomy_id;
	}
	public void setNcbi_taxonomy_id(String ncbi_taxonomy_id) {
		this.ncbi_taxonomy_id = ncbi_taxonomy_id;
	}
	public String getOrganism_common_name() {
		return organism_common_name;
	}
	public void setOrganism_common_name(String organism_common_name) {
		this.organism_common_name = organism_common_name;
	}
	public String getOrganism_scientific() {
		return organism_scientific;
	}
	public void setOrganism_scientific(String organism_scientific) {
		this.organism_scientific = organism_scientific;
	}
	public String getStrain() {
		return strain;
	}
	public void setStrain(String strain) {
		this.strain = strain;
	}
}
