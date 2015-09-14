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

/** A simple bean to store disulfide bridge information, the SSBOND records in the PDB files.
 *
 * The two residues specified here are CYS residues that form a Disulfide bridge.
 *
 *
 *
 * @author Andreas Prlic
 *
 */
public interface SSBond extends PDBRecord, Serializable, Cloneable {
    @Override
	public String toPDB();

    /**
     * append the PDB representation of this SSBOND to the provided StringBUffer
     *
     * @param buf a StringBuffer to print the PDB representation to
     */
    @Override
	public void toPDB(StringBuffer buf);

    public String getInsCode1();

    public void setInsCode1(String insCode1);

    public String getInsCode2();

    public void setInsCode2(String insCode2);

    /**
     * set serial number of this SSBOND in PDB file
     *
     * @return the serial number
     */
    public int getSerNum();

    /**
     * get serial number of this SSBOND in PDB file
     *
     * @param serNum
     */
    public void setSerNum(int serNum);

    public String getChainID1();

    public void setChainID1(String chainID1);

    public String getChainID2();

    public void setChainID2(String chainID2);

    /**
     * get residue number for first CYS. number and insertion code are joint
     * together.
     *
     * @return the residue number of the first CYS.
     *
     */
    public String getResnum1();

    public void setResnum1(String resnum1);

    /**
     * get residue number for second CYS. number and insertion code are joint
     * together.
     *
     * @return the residue number of the second CYS.
     *
     */
    public String getResnum2();

    public void setResnum2(String resnum2);
    
    public SSBond clone() ;
}
