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
package org.biojava.nbio.structure;

import junit.framework.TestCase;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;

public class TestDownloadChemCompProvider extends TestCase{
	
	public void testProtectedIDs(){
		
		String id = "CON";
		
		DownloadChemCompProvider prov = new DownloadChemCompProvider();
		ChemComp cc = prov.getChemComp(id);
		assertNotNull(cc);
		
		assertEquals(cc.getId(), id);
	}

}
