package demo;

import java.io.IOException;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.secstruc.SecStrucInfo;
import org.biojava.nbio.structure.secstruc.SecStrucPred;

/**
 * Demonstration of how to use the Secondary Structure Prediction (DSSP)
 * implementation in BioJava.
 * 
 * @author Aleix Lafita
 *
 */
public class DemoSecStrucPred {

    public static void main(String[] args) 
    		throws IOException, StructureException {
	
    	String pdbID = "5pti";
    	
        AtomCache cache = new AtomCache();
        
        //Load structure without any SS assignment
        Structure s = cache.getStructure(pdbID);
        
        //Predict and assign the SS of the Structure
        SecStrucPred ssp = new SecStrucPred();
        ssp.predict(s, true);
	
	    for (Chain c : s.getChains()) {
	        for (Group g: c.getAtomGroups()){
	
	            if (g.hasAminoAtoms()){
	
	                SecStrucInfo ss = 
	                		(SecStrucInfo) g.getProperty(Group.SEC_STRUC);
	
	                System.out.println(c.getChainID() + 
	                		" " + g.getResidueNumber() + " " 
	                		+ g.getPDBName() + " -> " + ss);
	            }
	        }
	    }
    }
}
