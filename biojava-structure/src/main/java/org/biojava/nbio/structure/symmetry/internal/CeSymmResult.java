package org.biojava.nbio.structure.symmetry.internal;

import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;

/**
 * This class stores all the relevant information of an internal symmetry
 * result obtained with CeSymm. The purpose is to carry all the information
 * packed during the calculations and return a single object.
 * 
 * @author Aleix Lafita
 *
 */
public class CeSymmResult {

	private MultipleAlignment multipleAlignment;
	private AFPChain selfAlignment;
	
	private CESymmParameters params = new CESymmParameters();
	private SymmetryAxes axes;
	private QuatSymmetryResults PG;
	
	private int symmOrder = 1;
	private boolean refined = false;
	private SymmetryType type;
	
	public MultipleAlignment getMultipleAlignment() {
		return multipleAlignment;
	}
	public void setMultipleAlignment(MultipleAlignment multipleAlignment) {
		this.multipleAlignment = multipleAlignment;
	}
	public AFPChain getSelfAlignment() {
		return selfAlignment;
	}
	public void setSelfAlignment(AFPChain selfAlignment) {
		this.selfAlignment = selfAlignment;
	}
	public CESymmParameters getParams() {
		return params;
	}
	public void setParams(CESymmParameters params) {
		this.params = params;
	}
	public SymmetryAxes getAxes() {
		return axes;
	}
	public void setAxes(SymmetryAxes axes) {
		this.axes = axes;
	}
	public int getSymmOrder() {
		return symmOrder;
	}
	public void setSymmOrder(int symmOrder) {
		this.symmOrder = symmOrder;
	}
	public boolean isRefined() {
		return refined;
	}
	public void setRefined(boolean refined) {
		this.refined = refined;
	}
	public QuatSymmetryResults getPG() {
		return PG;
	}
	public void setPG(QuatSymmetryResults pG) {
		PG = pG;
	}
	public SymmetryType getType() {
		return type;
	}
	public void setType(SymmetryType type) {
		this.type = type;
	}
	
}
