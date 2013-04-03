package org.biojava.bio.structure;

import org.biojava.bio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;

import junit.framework.TestCase;

public class TestDownloadChemCompProvider extends TestCase{
	
	public void testProtectedIDs(){
		
		String id = "CON";
		
		DownloadChemCompProvider prov = new DownloadChemCompProvider();
		ChemComp cc = prov.getChemComp(id);
		assertNotNull(cc);
		
		assertEquals(cc.getId(), id);
	}

}
