package org.biojava.bio.structure.io.mmcif.model;

/** A bean for the Pdbx_entity_nonpoly category.
 *
 * @author Andreas Prlic
 * @since 1.7
 */
public class PdbxEntityNonPoly {
	String entity_id;
	String name;
	String comp_id;
	public String getEntity_id() {
		return entity_id;
	}
	public void setEntity_id(String entity_id) {
		this.entity_id = entity_id;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public String getComp_id() {
		return comp_id;
	}
	public void setComp_id(String comp_id) {
		this.comp_id = comp_id;
	}

}
