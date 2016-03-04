package org.biojava.nbio.structure;

import static org.junit.Assert.*;

import java.io.IOException;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.quaternary.io.BioUnitDataProvider;
import org.biojava.nbio.structure.quaternary.io.BioUnitDataProviderFactory;
import org.biojava.nbio.structure.quaternary.io.MmCifBiolAssemblyProvider;
import org.junit.Test;

public class TestBioAssemblyIdentifier {


	@Test
	public void test() throws IOException, StructureException {
		// Save global state
		Class<? extends BioUnitDataProvider> prevProvider = BioUnitDataProviderFactory.getBioUnitDataProviderClass();
		BioUnitDataProviderFactory.setBioUnitDataProvider(MmCifBiolAssemblyProvider.class);
		
		AtomCache cache = new AtomCache();
		BioAssemblyIdentifier id;
		Structure s;
		
		// first assembly
		id = new BioAssemblyIdentifier("BIO:2ehz:1");
		s = cache.getStructure(id);
		assertEquals("Number of models",8, s.nrModels());
		assertEquals("Number of chains per model",1,s.getChains(0).size());
		// equivalent
		id = new BioAssemblyIdentifier("BIO:2ehz");
		s = cache.getStructure(id);
		assertEquals("Number of models",8, s.nrModels());
		assertEquals("Number of chains per model",1,s.getChains(0).size());
		// No second
		id = new BioAssemblyIdentifier("BIO:2ehz:2");
		try {
		s = cache.getStructure(id);
		fail("Expected exception for invalid assembly number");
		} catch( StructureException e) {}
		// AU
		id = new BioAssemblyIdentifier("BIO:2ehz:0");
		s = cache.getStructure(id);
		assertEquals("Number of models",1, s.nrModels());
		assertEquals("Number of chains per model",1,s.getChains(0).size());

		BioUnitDataProviderFactory.setBioUnitDataProvider(prevProvider);
	}

}
