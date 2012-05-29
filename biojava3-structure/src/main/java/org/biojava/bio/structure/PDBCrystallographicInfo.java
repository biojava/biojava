package org.biojava.bio.structure;

import java.io.Serializable;

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
	/**
	 * 
	 */
	private static final long serialVersionUID = -7949886749566087669L;
	private float a = 1.0f;
	private float b = 1.0f;
	private float c = 1.0f;
	private float alpha = 90.0f;
	private float beta = 90.0f;
	private float gamma = 90.0f;
	private String spaceGroup = "";
	private int z;
	
	/**
	 * @return the unit cell parameter a
	 */
	public float getA() {
		return a;
	}
	/**
	 * @param a the unit cell parameter a to set
	 */
	public void setA(float a) {
		this.a = a;
	}
	/**
	 * @return the unit cell parameter b
	 */
	public float getB() {
		return b;
	}
	/**
	 * @param b the unit cell parameter b to set
	 */
	public void setB(float b) {
		this.b = b;
	}
	/**
	 * @return the unit cell parameter c
	 */
	public float getC() {
		return c;
	}
	/**
	 * @param c the unit cell parameter c to set
	 */
	public void setC(float c) {
		this.c = c;
	}
	/**
	 * @return the unit cell parameter alpha (degrees)
	 */
	public float getAlpha() {
		return alpha;
	}
	/**
	 * @param alpha the unit cell parameter alpha to set
	 */
	public void setAlpha(float alpha) {
		this.alpha = alpha;
	}
	/**
	 * @return the unit cell parameter beta (degrees)
	 */
	public float getBeta() {
		return beta;
	}
	/**
	 * @param beta the unit cell parameter beta to set
	 */
	public void setBeta(float beta) {
		this.beta = beta;
	}
	/**
	 * @return the unit cell parameter gamma (degrees)
	 */
	public float getGamma() {
		return gamma;
	}
	/**
	 * @param gamma the unit cell parameter gamma to set
	 */
	public void setGamma(float gamma) {
		this.gamma = gamma;
	}
	/**
	 * @return the spaceGroup
	 */
	public String getSpaceGroup() {
		return spaceGroup;
	}
	/**
	 * @param spaceGroup the spaceGroup to set
	 */
	public void setSpaceGroup(String spaceGroup) {
		this.spaceGroup = spaceGroup;
	}
	/**
	 * @return the z
	 */
	public int getZ() {
		return z;
	}
	/**
	 * @param z the z to set
	 */
	public void setZ(int z) {
		this.z = z;
	}
	/**
	 * Returns true if structure was solved by a crystallographic method, i.e.,
	 * the values returned by this class are meaningful.
	 * @return the crystallographic
	 */
	public boolean isCrystallographic() {
		return spaceGroup.length() > 0;
	}

	
}
