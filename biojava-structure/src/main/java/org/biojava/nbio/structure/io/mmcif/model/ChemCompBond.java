/**
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
 * Created on Feb 5, 2013
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.nbio.structure.io.mmcif.model;

import java.io.Serializable;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/*
 * _chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
_chem_comp_bond.pdbx_ordinal
 */
public class ChemCompBond implements Serializable {

	private static final long serialVersionUID = 5905371029161975421L;

	private static final Logger logger = LoggerFactory.getLogger(ChemCompBond.class);

	String comp_id;
	String atom_id_1;
	String atom_id_2;
	String value_order;
	String pdbx_aromatic_flag;
	String pdbx_stereo_config;
	String pdbx_ordinal;
	public String getComp_id() {
		return comp_id;
	}
	public void setComp_id(String comp_id) {
		this.comp_id = comp_id;
	}
	public String getAtom_id_1() {
		return atom_id_1;
	}
	public void setAtom_id_1(String atom_id_1) {
		this.atom_id_1 = atom_id_1;
	}
	public String getAtom_id_2() {
		return atom_id_2;
	}
	public void setAtom_id_2(String atom_id_2) {
		this.atom_id_2 = atom_id_2;
	}
	public String getValue_order() {
		return value_order;
	}
	public void setValue_order(String value_order) {
		this.value_order = value_order;
	}
	public String getPdbx_aromatic_flag() {
		return pdbx_aromatic_flag;
	}
	public void setPdbx_aromatic_flag(String pdbx_aromatic_flag) {
		this.pdbx_aromatic_flag = pdbx_aromatic_flag;
	}
	public String getPdbx_stereo_config() {
		return pdbx_stereo_config;
	}
	public void setPdbx_stereo_config(String pdbx_stereo_config) {
		this.pdbx_stereo_config = pdbx_stereo_config;
	}
	public String getPdbx_ordinal() {
		return pdbx_ordinal;
	}
	public void setPdbx_ordinal(String pdbx_ordinal) {
		this.pdbx_ordinal = pdbx_ordinal;
	}

	/**
	 * Converts this ChemCompBond's value_order attribute into an int using the
	 * conversion:
	 *
	 * <pre>
	 * 	SING -> 1
	 * 	DOUB -> 2
	 * 	TRIP -> 3
	 * 	QUAD -> 4
	 * </pre>
	 *
	 * Any other values will return -1.
	 * <p>
	 * (Source:
	 * http://mmcif.rcsb.org/dictionaries/mmcif_mdb.dic/Items/_chem_comp_bond.
	 * value_order.html)
	 *
	 * @return the numerical value of this ChemCompBond's bond order, or -1 if
	 *         the value is non-numeric or unknown.
	 */
	public int getNumericalBondOrder() {
		if (value_order.equals("SING")) {
			return 1;
		} else if (value_order.equals("DOUB")) {
			return 2;
		} else if (value_order.equals("TRIP")) {
			return 3;
		} else if (value_order.equals("QUAD")) {
			return 4;
		} else {
			logger.error("Unknown or non-numeric value for value_order: "
					+ value_order);
			return -1;
		}
	}
}
