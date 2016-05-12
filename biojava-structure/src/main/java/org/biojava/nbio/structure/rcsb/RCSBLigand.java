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
 * Created on 2013-06-13
 * Created by Douglas Myers-Turnbull
 *
 * @since 3.0.6
 */
package org.biojava.nbio.structure.rcsb;

/**
 * Corresponds to a ligand in a {@code ligandInfo} XML file.
 *
 * @see <a href="http://www.pdb.org/pdb/software/rest.do#descPDB">RCSB RESTful</a>
 *
 * @author dmyerstu
 * @since 3.0.6
 */
public class RCSBLigand {

	private String formula;
	private String id;
	private String inChI;
	private String inChIKey;
	private String name;
	private String smiles;
	private String type;
	private Double weight;

	public String getFormula() {
		return formula;
	}

	public String getId() {
		return id;
	}

	public String getInChI() {
		return inChI;
	}

	public String getInChIKey() {
		return inChIKey;
	}

	public String getName() {
		return name;
	}

	public String getSmiles() {
		return smiles;
	}

	public String getType() {
		return type;
	}

	public Double getWeight() {
		return weight;
	}

	public void setFormula(String formula) {
		this.formula = formula;
	}

	public void setId(String id) {
		this.id = id;
	}

	public void setInChI(String inChI) {
		this.inChI = inChI;
	}

	public void setInChIKey(String inChIKey) {
		this.inChIKey = inChIKey;
	}

	public void setName(String name) {
		this.name = name;
	}

	public void setSmiles(String smiles) {
		this.smiles = smiles;
	}

	public void setType(String type) {
		this.type = type;
	}

	public void setWeight(Double weight) {
		this.weight = weight;
	}

}
