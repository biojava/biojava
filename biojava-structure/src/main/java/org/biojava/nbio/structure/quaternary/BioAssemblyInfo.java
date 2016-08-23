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
package org.biojava.nbio.structure.quaternary;

import java.io.Serializable;
import java.util.List;

/**
 * Representation of a Biological Assembly annotation as provided by the PDB.
 * Contains all the information required to build the Biological Assembly from
 * the asymmetric unit.
 * Note that the PDB allows for 1 or more Biological Assemblies for a given entry. They
 * are identified by the id field.
 *
 * @author Jose Duarte
 */
public class BioAssemblyInfo implements Serializable {


	private static final long serialVersionUID = 1L;

	private int id;
	private List<BiologicalAssemblyTransformation> transforms;
	private int macromolecularSize;

	/**
	 * Empty constructor
	 */
	public BioAssemblyInfo() {

	}

	/**
	 * The identifier for this Biological Assembly, from 1 to n
	 * @return
	 */
	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	/**
	 * Return the list of {@link BiologicalAssemblyTransformation}s needed to generate
	 * the biological assembly. There is one transformation per internal chain id.
	 * @return
	 */
	public List<BiologicalAssemblyTransformation> getTransforms() {
		return transforms;
	}

	public void setTransforms(List<BiologicalAssemblyTransformation> transforms) {
		this.transforms = transforms;
	}

	/**
	 * Returns the macromolecular size of this biological assembly, i.e.
	 * the number of polymeric chains (protein or nucleotide chains, not sugars) 
	 * in the biological assembly.
	 * @return
	 */
	public int getMacromolecularSize() {
		return macromolecularSize;
	}

	public void setMacromolecularSize(int macromolecularSize) {
		this.macromolecularSize = macromolecularSize;
	}

	@Override
	public String toString() {
		return "[BioAssembly "+id+": "+transforms.size()+" transforms, mm size "+macromolecularSize+"]";
	}
}
