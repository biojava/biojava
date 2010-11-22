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

import java.util.regex.Matcher;
import java.util.regex.Pattern;


/** This class uniquely describes an atom
 * 
 * @author Andreas Prlic
 *
 */
public class AtomInfo {

    String chainId;
    String atomName;
   
    String residueName;
    String residueNumber;
    int modelNumber;
    
    private static Pattern inscodePatter ;
	static {
		inscodePatter = Pattern.compile("([0-9]+)([a-zA-Z]*)?");
	}
    
    public AtomInfo() {
        super();

    }
    
    public static AtomInfo fromString(String atomInfo){
    	return AtomInfoParser.parse(atomInfo);
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

    /** Including insertion code
     * 
     * @param residueName
     */
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
		String aa3 = "";
		boolean printResName = true;

		String chain1 ="";
		String res1 = "";

		aa3 = residueName;				
		res1 = residueNumber;
		chain1 = chainId;
		
		StringBuffer buf = new StringBuffer();
		if ( printResName) {
			if ( !aa3.equals("")){
				buf.append("[");
				buf.append(aa3);
				buf.append("]");
			}
		}
		if ( ! res1.equals("")) {

			// let's check if there is an insertion code...
			Matcher matcher = inscodePatter.matcher(res1);

			boolean found = matcher.find();
			if ( ! found) {
				System.err.println("JmolTools: could not parse the residue number string " + res1);
				buf.append(res1);
			} else {
				String residueNumber = matcher.group(1);
				String insCode = matcher.group(2);
				buf.append(residueNumber);
				if ( insCode != null && ! ( insCode.equals(""))) {
					buf.append("^");
					buf.append(insCode);
				}								
			}

		}

		if ( ! chain1.equals("")){
			buf.append(":");
			buf.append(chain1);
		}
		
		if ( atomName != null) {
			buf.append(".");
			buf.append(atomName);
		}
		if ( modelNumber > 0) {
			buf.append("/");
			buf.append(modelNumber);
		}
		return buf.toString();
	}




}
