/*
 *                  BioJava development code
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
 * Created on Jul 25, 2006
 *
 */
package org.biojava.bio.structure.align.gui.jmol;

public class AtomInfo {

    String chainId;
    String atomName;
   
    String residueName;
    String residueNumber;
    int modelNumber;
    
    public AtomInfo() {
        super();

    }
    
    public int getModelNumber() {
        return modelNumber;
    }

    public void setModelNumber(int modelNumber) {
        this.modelNumber = modelNumber;
    }

    public String getResidueName() {
        return residueName;
    }

    public void setResidueName(String residueName) {
        this.residueName = residueName;
    }

    public String getResidueNumber() {
        return residueNumber;
    }

    public void setResidueNumber(String residueNumber) {
        this.residueNumber = residueNumber;
    }

    public String getChainId() {
        return chainId;
    }



    public void setChainId(String chainId) {
        this.chainId = chainId;
    }



    public String getAtomName() {
        return atomName;
    }



    public void setAtomName(String name) {
        this.atomName = name;
    }

	@Override
	public String toString() {
		return "AtomInfo [atomName=" + atomName + ", chainId=" + chainId
				+ ", modelNumber=" + modelNumber + ", residueName="
				+ residueName + ", residueNumber=" + residueNumber + "]";
	}




}
