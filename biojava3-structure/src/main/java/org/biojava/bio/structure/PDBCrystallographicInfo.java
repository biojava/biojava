package org.biojava.bio.structure;

import java.io.Serializable;

import org.biojava.bio.structure.xtal.CrystalCell;
import org.biojava.bio.structure.xtal.SpaceGroup;
import org.biojava.bio.structure.xtal.SymoplibParser;

/**
 * A class to hold crystallographic information about a PDB structure. The information
 * is only meaningful if the space group is defined. Use the method isCrystallographic()
 * to check if the associated PDB structure is of crystallographic origin. For example, these
 * values are meaningless for NMR structures.
 * 
 * @author Peter Rose
 *
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
	
	@Deprecated
	/**
	 * @param a the unit cell parameter a to set
	 */
	public void setA(float a) {
		this.cell.setA(a);
	}
	
	/**
	 * @return the unit cell parameter b
	 */
	public float getB() {
		return (float)cell.getB();
	}
	
	@Deprecated
	/**
	 * @param b the unit cell parameter b to set
	 */
	public void setB(double b) {
		cell.setB(b);
	}
	
	/**
	 * @return the unit cell parameter c
	 */
	public float getC() {
		return (float)cell.getC();
	}
	
	@Deprecated
	/**
	 * @param c the unit cell parameter c to set
	 */
	public void setC(float c) {
		cell.setC(c);
	}
	
	/**
	 * @return the unit cell parameter alpha (degrees)
	 */
	public float getAlpha() {
		return (float)cell.getAlpha();
	}
	
	@Deprecated
	/**
	 * @param alpha the unit cell parameter alpha to set
	 */
	public void setAlpha(float alpha) {
		cell.setAlpha(alpha);
	}
	
	/**
	 * @return the unit cell parameter beta (degrees)
	 */
	public float getBeta() {
		return (float)cell.getBeta();
	}
	
	@Deprecated
	/**
	 * @param beta the unit cell parameter beta to set
	 */
	public void setBeta(float beta) {
		cell.setBeta(beta);
	}
	
	/**
	 * @return the unit cell parameter gamma (degrees)
	 */
	public float getGamma() {
		return (float)cell.getGamma();
	}
	
	@Deprecated
	/**
	 * @param gamma the unit cell parameter gamma to set
	 */
	public void setGamma(float gamma) {
		cell.setGamma(gamma); 
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
	 * @return the spaceGroup
	 */
	public String getSpaceGroup() {
		if (sg==null) return null;
		return sg.getShortSymbol();
	}
	
	//public SpaceGroup getSpaceGroup() {
	//	return sg;
	//}
	
	@Deprecated
	/**
	 * Use #setSpaceGroup(SpaceGroup) instead
	 * @param spaceGroup the spaceGroup to set
	 */
	public void setSpaceGroup(String spaceGroup) {
		this.sg = SymoplibParser.getSpaceGroup(spaceGroup);
	}
	
	/**
	 * Set the space group
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
	
	@Deprecated
	/**
	 * Returns true if structure was solved by a crystallographic method, i.e.,
	 * the values returned by this class are meaningful.
	 * Use {@link Structure.isCrystallographic} instead
	 * @return 
	 */
	public boolean isCrystallographic() {
		return sg!=null;
	}

	
}
