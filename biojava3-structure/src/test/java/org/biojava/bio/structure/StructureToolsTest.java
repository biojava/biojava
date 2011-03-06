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
 * Created on Jun 8, 2007
 *
 */
package org.biojava.bio.structure;

import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileParser;

import junit.framework.TestCase;

public class StructureToolsTest extends TestCase {

    Structure structure, structure2, structure3;

    @Override
    protected void setUp() throws IOException
    {
    	InputStream inStream = this.getClass().getResourceAsStream("/5pti.pdb");
    	assertNotNull(inStream);


    	PDBFileParser pdbpars = new PDBFileParser();
    	FileParsingParameters params = new FileParsingParameters();
    	params.setAlignSeqRes(false);
    	pdbpars.setFileParsingParameters(params);

    	structure = pdbpars.parsePDBFile(inStream) ;

    	assertNotNull(structure);

    	assertEquals("structure does not contain one chain ", 1 ,structure.size());

    	inStream.close();

    	// Load structure2
    	inStream = this.getClass().getResourceAsStream("/1lnl.pdb");
    	assertNotNull(inStream);

    	structure2 = pdbpars.parsePDBFile(inStream) ;

    	assertNotNull(structure2);

    	assertEquals("structure does not contain 3 chains ", 3 ,structure2.size());
    	
    	inStream.close();

    	// Load structure2
    	inStream = this.getClass().getResourceAsStream("/1a4w.pdb");
    	assertNotNull(inStream);

    	structure3 = pdbpars.parsePDBFile(inStream) ;

    	assertNotNull(structure3);

    	assertEquals("structure does not contain 3 chains ", 3 ,structure3.size());
    	
    	inStream.close();
    }


    public void testGetCAAtoms(){
        Atom[] cas = StructureTools.getAtomCAArray(structure);
        assertEquals("did not find the expected number of Atoms (58), but got " + cas.length,58,cas.length);
    }

    public void testGetNrAtoms(){
        int length = StructureTools.getNrAtoms(structure);
        assertEquals("did not find the expected number of Atoms (1070), but got " + length,1070,length);


    }
    
    public void testGetSubRanges() throws StructureException {
    	String range;
    	Structure substr;
    	Chain chain;
    	
    	// normal substructures
    	range = "A:3-7";
    	substr = StructureTools.getSubRanges(structure2, range);
    	assertEquals("Wrong number of chains in "+range, 1, substr.size());

    	chain = substr.getChain(0);

    	assertEquals("Did not find the expected number of residues in "+range, 5, chain.getAtomLength() );
    	
    	// full chains
    	range = "A";
    	substr = StructureTools.getSubRanges(structure2, range);
    	assertEquals("Wrong number of chains in "+range, 1, substr.size());

    	chain = substr.getChain(0);

    	assertEquals("Did not find the expected number of residues in "+range, 411, chain.getAtomLength() );
    	//assertEquals("subrange doesn't equal original chain A.", structure2.getChainByPDB("A"), chain);
    	
    	// full chains
    	range = "A:";
    	substr = StructureTools.getSubRanges(structure2, range);
    	assertEquals("Wrong number of chains in "+range, 1, substr.size());

    	chain = substr.getChain(0);

    	assertEquals("Did not find the expected number of residues in "+range, 411, chain.getAtomLength() );
    	//assertEquals("subrange doesn't equal original chain A.", structure2.getChainByPDB("A"), chain);
    	
    	// combined ranges
    	range = "A:3-7,B:8-12";
    	substr = StructureTools.getSubRanges(structure2, range);
    	assertEquals("Wrong number of chains in "+range, 2, substr.size());

    	chain = substr.getChain(0);
    	assertEquals("Did not find the expected number of residues in first chain of "+range, 5, chain.getAtomLength() );
    	
    	chain = substr.getChain(1);
    	assertEquals("Did not find the expected number of residues in second chain of "+range, 5, chain.getAtomLength() );

    	// combined ranges
    	range = "A,B:8-12";
    	substr = StructureTools.getSubRanges(structure2, range);
    	assertEquals("Wrong number of chains in "+range, 2, substr.size());

    	chain = substr.getChain(0);
    	assertEquals("Did not find the expected number of residues in first chain of "+range, 411, chain.getAtomLength() );
    	
    	chain = substr.getChain(1);
    	assertEquals("Did not find the expected number of residues in second chain of "+range, 5, chain.getAtomLength() );
    	
    	// parentheses
    	range = "(A:3-7)";
    	substr = StructureTools.getSubRanges(structure2, range);
    	assertEquals("Wrong number of chains in "+range, 1, substr.size());

    	chain = substr.getChain(0);
    	assertEquals("Did not find the expected number of residues in "+range, 5, chain.getAtomLength() );

    }
    
    /**
     * Test some subranges that we used to have problems with
     * @throws StructureException
     */
    public void testGetSubRangesExtended() throws StructureException {
    	String range;
    	Structure substr;
    	Chain chain;
    	
    	// negative indices
    	range = "A:-3-7";
    	substr = StructureTools.getSubRanges(structure2, range);
    	assertEquals("Wrong number of chains in "+range, 1, substr.size());

    	chain = substr.getChain(0);

    	// Note residue 0 is missing from 1lnl
    	assertEquals("Did not find the expected number of residues in "+range, 10, chain.getAtomLength() );
    	
    	// double negative indices
    	range = "A:-3--1";
    	substr = StructureTools.getSubRanges(structure2, range);
    	assertEquals("Wrong number of chains in "+range, 1, substr.size());

    	chain = substr.getChain(0);

    	assertEquals("Did not find the expected number of residues in "+range, 3, chain.getAtomLength() );
    	
    	// mixed indices
    	range = "A:-3-+1";
    	substr = StructureTools.getSubRanges(structure2, range);
    	assertEquals("Wrong number of chains in "+range, 1, substr.size());

    	chain = substr.getChain(0);

    	assertEquals("Did not find the expected number of residues in "+range, 4, chain.getAtomLength() );
    	
    	// positive indices
    	range = "A:+1-6";
    	substr = StructureTools.getSubRanges(structure2, range);
    	assertEquals("Wrong number of chains in "+range, 1, substr.size());

    	chain = substr.getChain(0);

    	assertEquals("Did not find the expected number of residues in "+range, 6, chain.getAtomLength() );
    	
    	
    	// whitespace
    	range = "A:3-7, B:8-12";
    	substr = StructureTools.getSubRanges(structure2, range);
    	assertEquals("Wrong number of chains in "+range, 2, substr.size());

    	chain = substr.getChain(0);
    	assertEquals("Did not find the expected number of residues in first chain of "+range, 5, chain.getAtomLength() );
    	
    	chain = substr.getChain(1);
    	assertEquals("Did not find the expected number of residues in second chain of "+range, 5, chain.getAtomLength() );
    	    
    }

    /**
     * Test insertion codes
     * @throws StructureException
     */
    public void testGetSubRangesInsertionCodes() throws StructureException {
    	String range;
    	Structure substr;
    	Chain chain;
    	
    	// range including insertion
    	range = "H:35-37"; //includes 36A
    	substr = StructureTools.getSubRanges(structure3, range);
    	assertEquals("Wrong number of chains in "+range, 1, substr.size());

    	chain = substr.getChain(0);

    	assertEquals("Did not find the expected number of residues in "+range, 4, chain.getAtomLength() );
    	

    	// end with insertion
    	range = "H:35-36A";
    	substr = StructureTools.getSubRanges(structure3, range);
    	assertEquals("Wrong number of chains in "+range, 1, substr.size());

    	chain = substr.getChain(0);

    	assertEquals("Did not find the expected number of residues in "+range, 3, chain.getAtomLength() );
    	
    	// begin with insertion
    	range = "H:36A-38"; //includes 36A
    	substr = StructureTools.getSubRanges(structure3, range);
    	assertEquals("Wrong number of chains in "+range, 1, substr.size());

    	chain = substr.getChain(0);

    	assertEquals("Did not find the expected number of residues in "+range, 3, chain.getAtomLength() );
    	
    	// within insertion
    	range = "L:14-14K"; //includes 36A
    	substr = StructureTools.getSubRanges(structure3, range);
    	assertEquals("Wrong number of chains in "+range, 1, substr.size());

    	chain = substr.getChain(0);

    	assertEquals("Did not find the expected number of residues in "+range, 12, chain.getAtomLength() );
    	
    	// within insertion
    	range = "L:14C-14J"; //includes 36A
    	substr = StructureTools.getSubRanges(structure3, range);
    	assertEquals("Wrong number of chains in "+range, 1, substr.size());

    	chain = substr.getChain(0);

    	assertEquals("Did not find the expected number of residues in "+range, 8, chain.getAtomLength() );
    }
    
    public void testGroupsWithinShell() {
    	//TODO
    }

}
