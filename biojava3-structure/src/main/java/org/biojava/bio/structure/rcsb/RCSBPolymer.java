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
 * Created on 2012-11-20
 * Created by Douglas Myers-Turnbull
 *
 * @since 3.0.6
 */
package org.biojava.bio.structure.rcsb;

import java.util.ArrayList;
import java.util.List;

/**
 * Corresponds to a polymer in a {@code describeMol} XML file.
 * 
 * @see <a href="http://www.pdb.org/pdb/software/rest.do#descPDB">RCSB RESTful</a>
 * 
 * @author dmyerstu
 * @since 3.0.6
 */
public class RCSBPolymer {

	private List<Character> chains;

	private String description;

	private String enzClass;

	private Integer index;

	private Integer length;

	private RCSBMacromolecule molecule;

	private List<String> synonyms;

	private RCSBTaxonomy taxonomy;

	private String type;

	private Double weight;

	public RCSBPolymer() {
		chains = new ArrayList<Character>();
		synonyms = new ArrayList<String>();
	}

	public List<Character> getChains() {
		return chains;
	}

	public String getDescription() {
		return description;
	}

	public String getEnzClass() {
		return enzClass;
	}

	public Integer getIndex() {
		return index;
	}

	public Integer getLength() {
		return length;
	}

	public RCSBMacromolecule getMolecule() {
		return molecule;
	}

	public List<String> getSynonyms() {
		return synonyms;
	}

	public RCSBTaxonomy getTaxonomy() {
		return taxonomy;
	}

	public String getType() {
		return type;
	}

	public Double getWeight() {
		return weight;
	}

	void addChain(char chain) {
		chains.add(chain);
	}

	void addSynonym(String synonym) {
		synonyms.add(synonym);
	}

	void setDescription(String description) {
		this.description = description;
	}

	void setEnzClass(String enzClass) {
		this.enzClass = enzClass;
	}

	void setIndex(Integer index) {
		this.index = index;
	}

	void setLength(Integer length) {
		this.length = length;
	}

	void setMolecule(RCSBMacromolecule molecule) {
		this.molecule = molecule;
	}

	void setTaxonomy(RCSBTaxonomy taxonomy) {
		this.taxonomy = taxonomy;
	}

	void setType(String string) {
		type = string;
	}

	void setWeight(Double weight) {
		this.weight = weight;
	}

}
