package org.biojava.nbio.structure.symmetry.internal;

import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;

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
	
	private CESymmParameters params;
	private SymmetryAxes axes;
	
	private int symmOrder;
	private boolean refined;
	
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
	
}
