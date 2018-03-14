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
package org.biojava.nbio.structure.io;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.junit.Test;

/**
 * This class will test parsing short CONECT lines.
 * @since Nov 30, 2016
 * @author larsonm
 */
public class TestShortLines {
	
	@Test
	public void testConect() throws IOException {
		PDBFileParser pdbPars = new PDBFileParser();
		FileParsingParameters params = pdbPars.getFileParsingParameters();
		params.setCreateAtomBonds(true);
		
		// CONECTS will be deprecated, but will we create bonds?
		// Like the LINK records, should BioJava create BondImpl when params.setCreateAtomBonds(true)?
		
		StringBuilder sb = new StringBuilder();
		sb.append("HETATM 2398  P   FAD A 500       8.398  46.448  73.490  1.00 13.51           P \n");
		sb.append("HETATM 2399  PA  FAD A 500       6.089  45.580  75.235  1.00 15.88           P \n");
		sb.append("HETATM 2400  O1P FAD A 500       7.908  47.684  72.869  1.00  4.19           O \n");
		sb.append("CONECT 2400 2398\n");
		String shortLine = sb.toString();
		Structure s;
		// Parse short
		try(InputStream is = new ByteArrayInputStream(shortLine.getBytes())) {
			s = pdbPars.parsePDBFile(is);
		}
		
		// After 4.2, CONECTS are deprecated, but there is not yet an implementation
		// describing how CONECTS will be replaced - will Bonds be created?
		// assertEquals(1, s.getConnections().size());
		assertNotNull(s); 
	}
	
	@Test
	public void testLINK() throws IOException {
		Structure s;
		PDBFileParser pdbPars = new PDBFileParser();
		FileParsingParameters params = pdbPars.getFileParsingParameters();
		params.setCreateAtomBonds(true);
		
		StringBuilder sb = new StringBuilder();
		sb.append("ATOM   2412  C21 2EG A   7       0.888  44.973  72.238  1.00 29.17           C \n");
		sb.append("ATOM   2413  C22 2EG B  19       0.888  44.973  72.238  1.00 29.17           C \n");
		//sb.append("LINK         C21 2EG A   7                 C22 2EG B  19     1555   1555  1.56 ");
		sb.append("LINK         C21 2EG A   7                 C22 2EG B  19\n");
		String shortLine = sb.toString();
		
		// Parse short
		try(InputStream is = new ByteArrayInputStream(shortLine.getBytes())) {
			s = pdbPars.parsePDBFile(is);
		}
		
		// Should be a bond present in the Atoms.
		Chain c = s.getChainByIndex(0, 0);
		Group g = c.getAtomGroups().get(0);
		Atom a = g.getAtom(0);
		assertEquals(1, a.getBonds().size());
	}
}
