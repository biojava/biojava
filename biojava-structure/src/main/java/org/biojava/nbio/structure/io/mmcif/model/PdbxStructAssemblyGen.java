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

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import java.io.Serializable;

@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class PdbxStructAssemblyGen implements Serializable{
	/**
	 *
	 */
	private static final long serialVersionUID = 6739568389242514332L;
	String assembly_id;
	String oper_expression;
	String asym_id_list;


	public String getAssembly_id() {
		return assembly_id;
	}
	public void setAssembly_id(String assembly_id) {
		this.assembly_id = assembly_id;
	}
	public String getOper_expression() {
		return oper_expression;
	}
	public void setOper_expression(String oper_expression) {
		this.oper_expression = oper_expression;
	}
	public String getAsym_id_list() {
		return asym_id_list;
	}
	public void setAsym_id_list(String asym_id_list) {
		this.asym_id_list = asym_id_list;
	}
	@Override
	public String toString() {
		return "PdbxStructAssemblyGen [assembly_id=" + assembly_id
				+ ", oper_expression=" + oper_expression + ", asym_id_list="
				+ asym_id_list + "]";
	}





}
