package org.biojava.bio.structure;

import java.io.Serializable;

import javax.vecmath.Matrix4d;

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
	
	/**
	 * Some PDB files contain NCS operators necessary to create the full AU.
	 * Usually this happens for viral proteins.
	 * See http://www.wwpdb.org/documentation/format33/sect8.html#MTRIXn .
	 * Note that the "given" operators 
	 * (iGiven field =1 in PDB format, "given" string in _struct_ncs_oper.code in mmCIF format) 
	 * are not stored.
	 */
	private Matrix4d[] ncsOperators;
	
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
	 * @return 0
	 * @deprecated As of 4.0, use {@link SpaceGroup#getMultiplicity()} and {@link Structure#size()}
	 */
	public int getZ() {		
		return 0;
	}
	
	/**
	 * Gets all symmetry transformation operators corresponding to this structure's space group 
	 * (including the identity, at index 0) expressed in the orthonormal basis. Using PDB axes 
	 * convention (NCODE=1).
	 * @return an array of size {@link SpaceGroup#getNumOperators()}
	 */	
	public Matrix4d[] getTransformationsOrthonormal() {
		Matrix4d[] transfs = new Matrix4d[this.getSpaceGroup().getNumOperators()];
		transfs[0] = new Matrix4d(this.getSpaceGroup().getTransformation(0)); // no need to transform the identity
		for (int i=1;i<this.getSpaceGroup().getNumOperators();i++) {
			transfs[i] = this.cell.transfToOrthonormal(this.getSpaceGroup().getTransformation(i));
		}
		return transfs;
	}
	
	/**
	 * Get the NCS operators.
	 * Some PDB files contain NCS operators necessary to create the full AU.
	 * Usually this happens for viral proteins.
	 * See http://www.wwpdb.org/documentation/format33/sect8.html#MTRIXn .
	 * Note that the "given" operators 
	 * (iGiven field =1 in PDB format, "given" string in _struct_ncs_oper.code in mmCIF format) 
	 * are not stored. 
	 * @return the operators or null if no operators present
	 */
	public Matrix4d[] getNcsOperators() {
		return ncsOperators;
	}

	/**
	 * Set the NCS operators.
	 * Some PDB files contain NCS operators necessary to create the full AU.
	 * Usually this happens for viral proteins.
	 * See http://www.wwpdb.org/documentation/format33/sect8.html#MTRIXn .
	 * Note that the "given" operators 
	 * (iGiven field =1 in PDB format, "given" string in _struct_ncs_oper.code in mmCIF format) 
	 * are not stored.
	 * @param ncsOperators 
	 */
	public void setNcsOperators(Matrix4d[] ncsOperators) {
		this.ncsOperators = ncsOperators;
	}
	
	
	@Override
	public String toString() {
		return "["+ 
				(sg==null?"no SG":sg.getShortSymbol())+" - "+
				
				(cell==null?"no Cell":
				String.format("%.2f %.2f %.2f, %.2f %.2f %.2f",
						cell.getA(),cell.getB(),cell.getC(),cell.getAlpha(),cell.getBeta(),cell.getGamma()) )+
				(ncsOperators==null? "" : String.format(" - %d NCS operators",ncsOperators.length) )+
				"]";
	}

	
}
