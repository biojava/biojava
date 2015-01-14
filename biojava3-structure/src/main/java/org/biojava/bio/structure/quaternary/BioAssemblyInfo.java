package org.biojava.bio.structure.quaternary;

import java.util.List;

/**
 * Representation of a Biological Assembly annotation as provided by the PDB.
 * Contains all the information required to build the Biological Assembly from 
 * the asymmetric unit.
 * Note that the PDB allows for 1 or more Biological Assemblies for a given entry. They
 * are identified by the id field.
 * 
 * @author duarte_j
 */
public class BioAssemblyInfo {
	
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
