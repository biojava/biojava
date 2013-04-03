package org.biojava.bio.structure.redmine;


import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.NucleotideImpl;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;

import junit.framework.TestCase;

/** test for https://redmine.open-bio.org/issues/3282
 * 
 * @author Andreas Prlic
 *
 */
public class test1DARSeqAlign extends TestCase {

	public void test1DAR(){
		AtomCache cache = new AtomCache();
		FileParsingParameters orig = cache.getFileParsingParams();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		params.setLoadChemCompInfo(true);
		
		cache.setFileParsingParams(params);
		
		
		try {
			Structure struc = cache.getStructure("1DAR");
			//System.out.println(struc);
			Chain c = struc.getChainByPDB("A");
			//System.out.println(c.getSeqResGroups());
			
			Group g = c.getGroupByPDB(ResidueNumber.fromString("692"));
			//System.out.println(g);
			//System.out.println(FileConvert.toPDB(g.getAtom(0)));
			
			Group g3 = c.getGroupByPDB(ResidueNumber.fromString("689"));
			//System.out.println(g3);
			//System.out.println(FileConvert.toPDB(g3.getAtom(0)));
			
			assertTrue(! c.getSeqResGroups().contains(g));
			
			assertTrue( g instanceof NucleotideImpl);
			
			assertTrue(g.getType().equals("nucleotide"));
			
			assertTrue( g3.getPDBName().equals("LYS"));
			assertTrue( c.getSeqResGroups().contains(g3));
			
			assertTrue( g3 instanceof AminoAcid);
					
			
		}
		catch (Exception e){
			fail(e.getMessage());
		}
		
		
		
		
		cache.setFileParsingParams(orig);
	}
}
