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
package org.biojava.nbio.structure.xtal;

import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.contact.BoundingBox;

import java.util.List;

import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

/**
 * A class to contain the BoundingBoxes of all polymeric molecules in a full unit cell.
 *
 * @author Jose Duarte
 *
 */
public class UnitCellBoundingBox {

	/**
	 * An array with dimensions numOperatorsSg x numPolyChainsAu to contain all
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
	private int numPolyChainsAu;

	public UnitCellBoundingBox(int numOperatorsSg, int numPolyChainsAu) {
		this.numOperatorsSg = numOperatorsSg;
		this.numPolyChainsAu = numPolyChainsAu;
		this.chainBbs = new BoundingBox[numOperatorsSg][numPolyChainsAu];
		this.auBbs = new BoundingBox[numOperatorsSg];
	}

	public void setBbs(Structure structure, Matrix4d[] ops, boolean includeHetAtoms) {

		setBb(structure, includeHetAtoms, 0);
		for (int i=1;i<ops.length;i++) {
			Structure sym = structure.clone();
			Calc.transform(sym, ops[i]);
			setBb(sym, includeHetAtoms, i);
		}

	}

	private void setBb(Structure s, boolean includeHetAtoms, int i) {
		chainBbs[i] = new BoundingBox[numPolyChainsAu];
		List<Chain> polyChains = s.getPolyChains();
		int j = 0;
		for (Chain polyChain : polyChains) {
			chainBbs[i][j] = new BoundingBox(StructureTools.getAllNonHCoordsArray(polyChain, includeHetAtoms));
			j++;
		}
		auBbs[i] = new BoundingBox(chainBbs[i]);
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
		UnitCellBoundingBox translatedBbs = new UnitCellBoundingBox(numOperatorsSg, numPolyChainsAu);

		for (int i=0; i<numOperatorsSg; i++) {
			for (int j = 0;j<numPolyChainsAu; j++) {
				translatedBbs.chainBbs[i][j] = new BoundingBox(this.chainBbs[i][j]);
				translatedBbs.chainBbs[i][j].translate(translation);
			}
			translatedBbs.auBbs[i] = new BoundingBox(translatedBbs.chainBbs[i]);
		}

		return translatedBbs;
	}

}
