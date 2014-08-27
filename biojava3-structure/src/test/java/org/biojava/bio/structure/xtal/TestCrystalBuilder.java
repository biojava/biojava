package org.biojava.bio.structure.xtal;

import static org.junit.Assert.*;

import java.io.IOException;

import org.biojava.bio.structure.ChainInterfaceList;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.structure.StructureIO;
import org.junit.Test;

public class TestCrystalBuilder {

	@Test
	public void test1NMR() throws IOException, StructureException { 

		// a monomer NMR entry: must have no interfaces
		
		AtomCache cache = new AtomCache();
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);
		Structure s1 = StructureIO.getStructure("1NMR");

		CrystalBuilder cb = new CrystalBuilder(s1);
		ChainInterfaceList interfaces = cb.getUniqueInterfaces(5.5);
		assertTrue(interfaces.size()==0);

	}
	
	@Test
	public void test1B8G() throws IOException, StructureException { 

		// a crystallographic entry: several interfaces
		
		AtomCache cache = new AtomCache();
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);
		Structure s1 = StructureIO.getStructure("1B8G");
		CrystalBuilder cb = new CrystalBuilder(s1);
		ChainInterfaceList interfaces = cb.getUniqueInterfaces(5.5);
		assertTrue(interfaces.size()>1);
		
		
	}
	
	@Test
	public void test2MFZ() throws IOException, StructureException { 

		// a dimer NMR entry: must have 1 interface
		
		AtomCache cache = new AtomCache();
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);
		Structure s1 = StructureIO.getStructure("2MFZ");
		CrystalBuilder cb = new CrystalBuilder(s1);
		ChainInterfaceList interfaces = cb.getUniqueInterfaces(5.5);
		assertTrue(interfaces.size()==1);
				
	}

	@Test
	public void test4MF8() throws IOException, StructureException { 

		// a crystallographic structure with protein+DNA: has only 3 prot-prot interfaces the rest are DNA-involving ones
		
		AtomCache cache = new AtomCache();
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);
		Structure s1 = StructureIO.getStructure("4MF8");
		CrystalBuilder cb = new CrystalBuilder(s1);
		ChainInterfaceList interfaces = cb.getUniqueInterfaces(5.5);
		assertTrue(interfaces.size()>3);
				
	}
}
