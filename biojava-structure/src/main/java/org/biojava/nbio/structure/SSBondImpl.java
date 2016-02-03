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
 * Created on Nov 14, 2007
 *
 */

package org.biojava.nbio.structure;

import java.io.Serializable;

/** 
 * A simple bean to store disulfide bridge information, the SSBOND records in the PDB files.
 *
 * The two residues specified here are CYS residues that form a Disulfide bridge.
 *
 *
 *
 * @author Andreas Prlic
 *
 */
public class SSBondImpl implements PDBRecord, Serializable, Cloneable, SSBond {

	private static final long serialVersionUID = -8663681100691188647L;

	private int serNum;
	
	private String chainID1;
	private String chainID2;
	private String resnum1;
	private String resnum2;
	private String insCode1;
	private String insCode2;

	public SSBondImpl(){
		serNum = 0;
	}

	@Override
	public String toPDB(){


		StringBuffer buf = new StringBuffer();
		toPDB(buf);
		return buf.toString();
	}

	/** 
	 * Append the PDB representation of this SSBOND to the provided StringBuffer
	 *
	 * @param buf a StringBuffer to print the PDB representation to
	 */
	@Override
	public void toPDB(StringBuffer buf){

		/*12 - 14        LString(3)      "CYS"        Residue name.
		16             Character       chainID1     Chain identifier.
		18 - 21        Integer         seqNum1      Residue sequence number.
		22             AChar           icode1       Insertion code.
		26 - 28        LString(3)      "CYS"        Residue name.
		30             Character       chainID2     Chain identifier.
		32 - 35        Integer         seqNum2      Residue sequence number.
		36             AChar           icode2       Insertion code.
		60 - 65        SymOP           sym1         Symmetry oper for 1st resid
		67 - 72        SymOP           sym2         Symmetry oper for 2nd resid
		*/
		//01234567890123456789012345678901234567890123456789012345678901234567890123456789
		//SSBOND   1 CYS      5    CYS     55                                     5PTI  67
		//SSBOND   2 CYS     14    CYS     38
		//SSBOND   3 CYS     30    CYS     51


		buf.append("SSBOND ");
		buf.append(String.format("%3d", serNum));
		buf.append(String.format(" CYS %s %4s%1s  ",chainID1,resnum1,insCode1));
		buf.append(String.format(" CYS %s %4s%1s  ",chainID2,resnum2,insCode2));
	}


        @Override
	public String getInsCode1() {
		return insCode1;
	}

        @Override
	public void setInsCode1(String insCode1) {
		this.insCode1 = insCode1;
	}

        @Override
	public String getInsCode2() {
		return insCode2;
	}

        @Override
	public void setInsCode2(String insCode2) {
		this.insCode2 = insCode2;
	}

	
        @Override
	public int getSerNum() {
		return serNum;
	}

	/** get serial number of this SSBOND in PDB file
	 *
	 * @param serNum
	 */
        @Override
	public void setSerNum(int serNum) {
		this.serNum = serNum;
	}

	@Override
	public SSBondImpl clone() {
		SSBondImpl nbond = new SSBondImpl();
		nbond.setChainID1(chainID1);
		nbond.setChainID2(chainID2);
		nbond.setResnum1(resnum1);
		nbond.setResnum2(resnum2);
		return nbond;
	}

        @Override
	public String getChainID1() {
		return chainID1;
	}
        @Override
	public void setChainID1(String chainID1) {
		this.chainID1 = chainID1;
	}
        @Override
	public String getChainID2() {
		return chainID2;
	}
        @Override
	public void setChainID2(String chainID2) {
		this.chainID2 = chainID2;
	}

	
        @Override
	public String getResnum1() {
		return resnum1;
	}
        @Override
	public void setResnum1(String resnum1) {
		this.resnum1 = resnum1;
	}

	
        @Override
	public String getResnum2() {
		return resnum2;
	}
        @Override
	public void setResnum2(String resnum2) {
		this.resnum2 = resnum2;
	}

	@Override
	public String toString() {
		String s = "[SSBOND:\n";

		s += "Atom 1:\n";
		s += "\tChain: " + chainID1 + "\n";
		s += "\tResidue #: " + resnum1 + "\n";
		s += "\tIns. Code: " + insCode1 + "\n";

		s += "Atom 2:\n";
		s += "\tChain: " + chainID2 + "\n";
		s += "\tResidue #: " + resnum2 + "\n";
		s += "\tIns. Code: " + insCode2 + "\n";

		s += "]";

		return s;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((chainID1 == null) ? 0 : chainID1.hashCode());
		result = prime * result + ((chainID2 == null) ? 0 : chainID2.hashCode());
		result = prime * result + ((insCode1 == null) ? 0 : insCode1.hashCode());
		result = prime * result + ((insCode2 == null) ? 0 : insCode2.hashCode());
		result = prime * result + ((resnum1 == null) ? 0 : resnum1.hashCode());
		result = prime * result + ((resnum2 == null) ? 0 : resnum2.hashCode());
		return result;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		SSBondImpl other = (SSBondImpl) obj;
		if (chainID1 == null) {
			if (other.chainID1 != null)
				return false;
		} else if (!chainID1.equals(other.chainID1))
			return false;
		if (chainID2 == null) {
			if (other.chainID2 != null)
				return false;
		} else if (!chainID2.equals(other.chainID2))
			return false;
		if (insCode1 == null) {
			if (other.insCode1 != null)
				return false;
		} else if (!insCode1.equals(other.insCode1))
			return false;
		if (insCode2 == null) {
			if (other.insCode2 != null)
				return false;
		} else if (!insCode2.equals(other.insCode2))
			return false;
		if (resnum1 == null) {
			if (other.resnum1 != null)
				return false;
		} else if (!resnum1.equals(other.resnum1))
			return false;
		if (resnum2 == null) {
			if (other.resnum2 != null)
				return false;
		} else if (!resnum2.equals(other.resnum2))
			return false;
		return true;
	}
}
