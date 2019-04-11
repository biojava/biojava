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
 * Created on Mar 29, 2014
 * Author: andreas
 *
 */

package org.biojava.nbio.structure.test.scop;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.junit.Assert;
import org.junit.Test;

import java.net.URL;

public class TestParsingScopDomainPDB {


	@Test
	public void testd1gena(){

		try {
			URL u = new URL("http://scop.berkeley.edu/astral/pdbstyle/ver=2.03&id=d1gena_&output=txt.");

			PDBFileParser p = new PDBFileParser();

			Structure s = p.parsePDBFile(u.openStream());

			//System.out.println(s);
			Assert.assertTrue(StructureTools.getAllAtomArray(s).length > 100);
		} catch (Exception e){

			e.printStackTrace();
			Assert.fail(e.getMessage());
		}
	}
}
