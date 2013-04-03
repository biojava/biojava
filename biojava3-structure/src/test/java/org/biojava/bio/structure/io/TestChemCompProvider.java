/**
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
 * Created on Mar 29, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure.io;

import java.io.IOException;

import junit.framework.TestCase;

import org.biojava.bio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.bio.structure.io.mmcif.ChemCompProvider;

/** Test case for https://redmine.open-bio.org/issues/3334
 * 
 * @author Andreas Prlic
 *
 */
public class TestChemCompProvider extends TestCase {

	public  void testChemCompProvider(){
		
		String pdbId = "1znf";
		
		FileParsingParameters params = new FileParsingParameters();
		params.setLoadChemCompInfo(true);
		
		PDBFileReader r = new PDBFileReader();
		r.setAutoFetch(true);
		r.setFileParsingParameters(params);
		
		ChemCompProvider prov = ChemCompGroupFactory.getChemCompProvider();
		
		//System.out.println(prov.getClass().getName());
		
		try {
			r.getStructureById(pdbId);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			fail(e.getMessage());
		}
		
		ChemCompProvider prov2 = ChemCompGroupFactory.getChemCompProvider();
		String name1 = prov.getClass().getName();
		String name2 = prov2.getClass().getName();
		assertEquals("The ChemCompProvider got modified from " + name1 + " to " + name2, name1, name2);
	
		
		params.setLoadChemCompInfo(false);
		r.setFileParsingParameters(params);
		
		try {
			r.getStructureById(pdbId);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			fail(e.getMessage());
		}
		
		ChemCompProvider prov3 = ChemCompGroupFactory.getChemCompProvider();
		String name3 = prov3.getClass().getName();
		
		assertEquals( "org.biojava.bio.structure.io.mmcif.ReducedChemCompProvider",name3);
	}
}
