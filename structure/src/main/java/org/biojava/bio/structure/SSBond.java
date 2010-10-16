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

package org.biojava.bio.structure;

import java.io.Serializable;

/** A simple bean to store disulfid bridge information, the SSBOND records in the PDB files.
 *
 * The two residues specified here are CYS residues that form a Disulfid bridge.
 *
 *
 *
 * @author Andreas Prlic
 *
 */
public class SSBond implements PDBRecord, Serializable{

	/**
    *
    */
   private static final long serialVersionUID = -8663681100691188647L;
   int serNum;
	String chainID1;
	String chainID2;
	String resnum1;
	String resnum2;
	String insCode1;
	String insCode2;

	public SSBond(){
		serNum = 0;
	}

	public String toPDB(){


		StringBuffer buf = new StringBuffer();
		toPDB(buf);
		return buf.toString();
	}

	/** append the PDB representation of this SSBOND to the provided StringBUffer
	 *
	 * @param buf a StringBuffer to print the PDB representation to
	 */
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


	public String getInsCode1() {
		return insCode1;
	}

	public void setInsCode1(String insCode1) {
		this.insCode1 = insCode1;
	}

	public String getInsCode2() {
		return insCode2;
	}

	public void setInsCode2(String insCode2) {
		this.insCode2 = insCode2;
	}

	/** set serial number of this SSBOND in PDB file
	 *
	 * @return the serial number
	 */
	public int getSerNum() {
		return serNum;
	}

	/** get serial number of this SSBOND in PDB file
	 *
	 * @param serNum
	 */
	public void setSerNum(int serNum) {
		this.serNum = serNum;
	}

	public SSBond clone() {
		SSBond nbond = new SSBond();
		nbond.setChainID1(chainID1);
		nbond.setChainID2(chainID2);
		nbond.setResnum1(resnum1);
		nbond.setResnum2(resnum2);
		return nbond;
	}

	public String getChainID1() {
		return chainID1;
	}
	public void setChainID1(String chainID1) {
		this.chainID1 = chainID1;
	}
	public String getChainID2() {
		return chainID2;
	}
	public void setChainID2(String chainID2) {
		this.chainID2 = chainID2;
	}

	/** get residue number for first CYS.
	 *  number and insertion code are joint together.
	 *
	 * @return the residue number of the first CYS.
	 *
	 */
	public String getResnum1() {
		return resnum1;
	}
	public void setResnum1(String resnum1) {
		this.resnum1 = resnum1;
	}

	/** get residue number for second CYS.
	 *  number and insertion code are joint together.
	 *
	 * @return the residue number of the second CYS.
	 *
	 */
	public String getResnum2() {
		return resnum2;
	}
	public void setResnum2(String resnum2) {
		this.resnum2 = resnum2;
	}



}
