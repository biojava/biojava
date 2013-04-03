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
 * created at Mar 4, 2008
 */
package org.biojava.bio.structure.io.mmcif.model;

/** A simple class to represent Entity records in mmCif files
 * 
 * @author Andreas Prlic
 *
 */
public class Entity {
	String id;
	
	String type; 
	String src_method; 
	String pdbx_description; 
	String formula_weight;
	String pdbx_number_of_molecules; 
	String details;
	String pdbx_mutation; 
	String pdbx_fragment;
	String pdbx_ec;

	public String toString(){
		StringBuffer buf = new StringBuffer();

        buf.append("Entity - id:").append(id);

        buf.append(" type:").append(type);
        buf.append(" src_method:").append(src_method);
        buf.append(" pdbx_description:").append(pdbx_description);
        buf.append(" formula_weight:").append(formula_weight);
        buf.append(" pdbx_number_f_molecules:").append(pdbx_number_of_molecules);
        buf.append(" details:").append(details);
        buf.append(" pdbx_mutation:").append(pdbx_mutation);
        buf.append(" pdbx_fragment:").append(pdbx_fragment);
        buf.append(" pdbx_ec:").append(pdbx_ec);
		
		return buf.toString();
	}
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	
	public String getType() {
		return type;
	}
	public void setType(String type) {
		this.type = type;
	}
	public String getSrc_method() {
		return src_method;
	}
	public void setSrc_method(String src_method) {
		this.src_method = src_method;
	}
	public String getPdbx_description() {
		return pdbx_description;
	}
	public void setPdbx_description(String pdbx_description) {
		this.pdbx_description = pdbx_description;
	}
	public String getFormula_weight() {
		return formula_weight;
	}
	public void setFormula_weight(String formula_weight) {
		this.formula_weight = formula_weight;
	}
	public String getPdbx_number_of_molecules() {
		return pdbx_number_of_molecules;
	}
	public void setPdbx_number_of_molecules(String pdbx_number_of_molecules) {
		this.pdbx_number_of_molecules = pdbx_number_of_molecules;
	}
	public String getDetails() {
		return details;
	}
	public void setDetails(String details) {
		this.details = details;
	}
	public String getPdbx_mutation() {
		return pdbx_mutation;
	}
	public void setPdbx_mutation(String pdbx_mutation) {
		this.pdbx_mutation = pdbx_mutation;
	}
	public String getPdbx_fragment() {
		return pdbx_fragment;
	}
	public void setPdbx_fragment(String pdbx_fragment) {
		this.pdbx_fragment = pdbx_fragment;
	}
	public String getPdbx_ec() {
		return pdbx_ec;
	}
	public void setPdbx_ec(String pdbx_ec) {
		this.pdbx_ec = pdbx_ec;
	}
	
	
}
