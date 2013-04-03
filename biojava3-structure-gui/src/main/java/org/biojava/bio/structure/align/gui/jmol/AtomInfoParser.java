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

import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.structure.align.gui.jmol.AtomInfo;

public class AtomInfoParser {

    static Logger    logger      = Logger.getLogger("org.biojava.spice");
    
    static Pattern pattern;
    
    static {
        
        String numberPattern = "\\[(.*)\\]([0-9^a-zA-Z]+)(:[a-zA-Z]*)?\\.([a-zA-Z]+)(/[0-9]*)?";
        pattern = Pattern.compile(numberPattern);
        
    }
    
    public AtomInfoParser() {
        super();

        
        
    }
    
    public static void main(String[]args){
    	String s1 = "[GLY]371:A.CA #2811";
    	String s2 = "[ASP]1^A:A.CA/2 #2";
        System.out.println(s1 + " got: " + AtomInfoParser.parse(s1));
        System.out.println(s2 + " got: " + AtomInfoParser.parse(s2));
    }
    
    
    /** parses e.g. 
     *  [MET]361:A.CA/1 #2843
     *  [GLY]339:A.CA #2573
     *  [ASN]44.CA #704
     *  
     * @param jmolAtomInfo
     * @return an AtomInfo
     */
    public static AtomInfo parse(String jmolAtomInfo){
       
        
        Matcher matcher = pattern.matcher(jmolAtomInfo);
           
        boolean found = matcher.find();
        if ( ! found) {
            logger.info("AtomInfoParser: could not parse the atomInfo string " + jmolAtomInfo);
            return new AtomInfo();
        }
        String residueName   = matcher.group(1);        
        String residueNumber = matcher.group(2);
        String chainId       = matcher.group(3);
        String atomName      = matcher.group(4);
        String modelNumber   = matcher.group(5);
        
   
        
        //System.out.println(jmolAtomInfo +" | " +  residueName + " number:" + residueNumber + " chain:" + chainId + 
        //		" atomName:" + atomName + " modelNumber:" + modelNumber );
        
        AtomInfo info = new AtomInfo();
        
        info.setAtomName(atomName);
        info.setResidueName(residueName);
        info.setResidueNumber(residueNumber.replaceAll("\\^",""));
        
        String ci = " ";
        if (chainId != null)
            ci = chainId.substring(1,chainId.length());        
        info.setChainId(ci);        
        
        int mn = 1;
        if ( modelNumber != null)
            mn = Integer.parseInt(modelNumber.substring(1,modelNumber.length()));            
        info.setModelNumber(mn);
       
        
        return info;
    }

}
