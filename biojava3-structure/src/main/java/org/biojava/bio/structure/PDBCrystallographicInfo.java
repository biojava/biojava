package org.biojava.bio.structure;

import java.io.Serializable;

import org.biojava.bio.structure.xtal.CrystalCell;
import org.biojava.bio.structure.xtal.SpaceGroup;

/**
 * A class to hold crystallographic information about a PDB structure.
 * 
 * @author Peter Rose
 * @author duarte_j
 */
public class PDBCrystallographicInfo implements Serializable {

	private static final long serialVersionUID = -7949886749566087669L;
	
	private CrystalCell cell;
	private SpaceGroup sg;
	
	private int z;
	
	public PDBCrystallographicInfo() {
		
	}
	
	/**
	 * @return the unit cell parameter a
	 */
	public float getA() {
		return (float)cell.getA();
	}
	
	/**
	 * @return the unit cell parameter b
	 */
	public float getB() {
		return (float)cell.getB();
	}
	
	/**
	 * @return the unit cell parameter c
	 */
	public float getC() {
		return (float)cell.getC();
	}
	
	/**
	 * @return the unit cell parameter alpha (degrees)
	 */
	public float getAlpha() {
		return (float)cell.getAlpha();
	}
	
	/**
	 * @return the unit cell parameter beta (degrees)
	 */
	public float getBeta() {
		return (float)cell.getBeta();
	}
	
	/**
	 * @return the unit cell parameter gamma (degrees)
	 */
	public float getGamma() {
		return (float)cell.getGamma();
	}
	
	/**
	 * Return the crystal cell
	 * @return
	 */
	public CrystalCell getCrystalCell() {
		return cell;
	}
	
	/**
	 * Set the crystal cell
	 * @param cell
	 */
	public void setCrystalCell(CrystalCell cell) {
		this.cell = cell;
	}
	
	/**
	 * Get the SpaceGroup
	 * @return the spaceGroup
	 */
	public SpaceGroup getSpaceGroup() {
		return sg;
	}
	
	/**
	 * Set the SpaceGroup
	 * @param spaceGroup
	 */
	public void setSpaceGroup(SpaceGroup spaceGroup) {
		this.sg = spaceGroup;
	}
	
	/**
	 * Return the z, i.e. the multiplicity of the space group times the number of chains in asymmetric unit
	 * @return the z
	 */
	public int getZ() {		
		return z;
	}
	
	/**
	 * Set the z 
	 * @param z
	 */
	public void setZ(int z) {
		this.z = z;
	}
	
	
}
