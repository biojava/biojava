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
package demo;

import java.io.IOException;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.secstruc.DSSPParser;
import org.biojava.nbio.structure.secstruc.SecStrucInfo;

/**
 * Demonstration of how to load a Structure with the SS information,
 * either from the PDB file annotation (Author's assignment) or from
 * the DSSP file in the PDB servers (DSSP assignment).
 * 
 * @author Aleix Lafita
 *
 */
public class DemoLoadSecStruc {
	
    public static void main(String[] args) 
    		throws IOException, StructureException {
    	
    	String pdbID = "5pti";
    	
        FileParsingParameters params = new FileParsingParameters();
        //Only change needed to the normal Structure loading
        params.setParseSecStruc(true); //this is false as DEFAULT

        AtomCache cache = new AtomCache();
        cache.setFileParsingParams(params);
        cache.setUseMmCif(false);

        //The loaded Structure contains the SS assigned by Author
        Structure s = cache.getStructure(pdbID);
        
        //If the more detailed DSSP prediction is required call this
        DSSPParser.fetch(pdbID, s, true);

        //Print the assignment residue by residue
        System.out.println("Residue assignment: ");
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