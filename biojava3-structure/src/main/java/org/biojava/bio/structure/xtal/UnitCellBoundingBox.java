package org.biojava.bio.structure.xtal;

import javax.vecmath.Vector3d;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.contact.BoundingBox;

/**
 * A class to contain the BoundingBoxes of all molecules in a full unit cell
 * 
 * @author duarte_j
 *
 */
public class UnitCellBoundingBox {

	/**
	 * An array with dimensions numOperatorsSg x numChainsAu to contain all 
	 * bounding boxes of all chains of all AUs in unit cell
	 * e.g. chainBbs[0] would be the bounding boxes for all chains in the original AU
	 */
	private BoundingBox[][] chainBbs;
	
	/**
	 * An array with dimensions numOperatorsSg to contain all bounding boxes of
	 * all AUs in unit cell
	 */
	private BoundingBox[] auBbs;
	
	private int numOperatorsSg; // i.e. multiplicity of space group
	private int numChainsAu;
	
	public UnitCellBoundingBox(int numOperatorsSg, int numChainsAu) {
		this.numOperatorsSg = numOperatorsSg;
		this.numChainsAu = numChainsAu;
		this.chainBbs = new BoundingBox[numOperatorsSg][numChainsAu];
		this.auBbs = new BoundingBox[numOperatorsSg];
	}
	
	/**
	 * Calculate the bounding boxes of all chains in original asymmetric unit.
	 * Considers all non-Hydrogen atoms of all residues
	 * @param s
	 * @param includeHetAtoms
	 */
	public void setOriginalAuBbs (Structure s, boolean includeHetAtoms) {
		
		chainBbs[0] = new BoundingBox[numChainsAu];
		
		int i = 0;
		for (Chain chain:s.getChains()) {
			chainBbs[0][i] = new BoundingBox(StructureTools.getAllNonHAtomArray(chain, includeHetAtoms));;
			i++;
		}
		
		auBbs[0] = new BoundingBox(chainBbs[0]);
	}
	
	/**
	 * Calculate the bounding boxes of all chains in the unit cell except for the original
	 * asymmetric unit (done with {@link #setOriginalAuBbs(Structure, boolean)}
	 * @param cell
	 * @param includeHetAtoms
	 */
	public void setAllNonAuBbs(Structure[] cell, boolean includeHetAtoms) {
		
		// the original AU (i=0) was already set by setOriginalAuBbs
		for (int i=1; i<numOperatorsSg; i++) {
			for (int j = 0;j<numChainsAu; j++) {
				chainBbs[i][j] = new BoundingBox(StructureTools.getAllNonHAtomArray(cell[i].getChain(j), includeHetAtoms));
			}
			auBbs[i] = new BoundingBox(chainBbs[i]);
		}	
	}
		
	/**
	 * Get the chain BoundingBox for the given cell index (cellIdx=0 would be original AU)
	 * and chain index
	 * @param cellIdx
	 * @param chainIdx
	 * @return
	 */
	public BoundingBox getChainBoundingBox(int cellIdx, int chainIdx) {
		return chainBbs[cellIdx][chainIdx];
	}
	
	/**
	 * Get the AU BoundingBox for the given cell index (cellIdx=0 would be original AU)
	 * The AU BoundingBox is the BoundingBox that bounds all chains belonging to the AU
	 * @param cellIdx
	 * @return
	 */
	public BoundingBox getAuBoundingBox(int cellIdx) {		
		return auBbs[cellIdx];		
	}
	
	/**
	 * Returns a new BoundingBoxes object containing the same bounds as this 
	 * BoundingBoxes object translated by the given translation
	 * @param translation
	 * @return
	 */
	public UnitCellBoundingBox getTranslatedBbs(Vector3d translation) {
		UnitCellBoundingBox translatedBbs = new UnitCellBoundingBox(numOperatorsSg, numChainsAu);
		
		for (int i=0; i<numOperatorsSg; i++) {
			for (int j = 0;j<numChainsAu; j++) {
				translatedBbs.chainBbs[i][j] = new BoundingBox(this.chainBbs[i][j]);
				translatedBbs.chainBbs[i][j].translate(translation);
			}
			translatedBbs.auBbs[i] = new BoundingBox(translatedBbs.chainBbs[i]);
		}
		
		return translatedBbs;
	}
}
