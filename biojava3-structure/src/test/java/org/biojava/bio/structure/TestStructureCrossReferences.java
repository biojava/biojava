package org.biojava.bio.structure;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.contact.StructureInterface;
import org.biojava.bio.structure.contact.StructureInterfaceList;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.xtal.CrystalBuilder;
import org.biojava3.structure.StructureIO;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class TestStructureCrossReferences {

	private static final Logger logger = LoggerFactory.getLogger(TestStructureCrossReferences.class);
	
			
	@Test
	public void testCrossReferencesMmCif() throws IOException, StructureException {
		
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(false); 
		cache.setFileParsingParams(params);
		
		StructureIO.setAtomCache(cache); 
		
		Structure structure = StructureIO.getStructure("1smt");
		
		System.out.println("Testing references in mmCIF loading with NO alignSeqRes");
		doFullTest(structure);
		
	}
	
	@Test
	public void testCrossReferencesMmCifAlignSeqRes() throws IOException, StructureException {
		
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true); 
		cache.setFileParsingParams(params);
		
		StructureIO.setAtomCache(cache); 
		
		Structure structure = StructureIO.getStructure("1smt");
		
		System.out.println("Testing references in mmCIF loading with alignSeqRes");
		doFullTest(structure);
		
	}

	
	@Test
	public void testCrossReferencesPdb() throws IOException, StructureException {
		
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(false);
		
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(false); 
		cache.setFileParsingParams(params);
		
		StructureIO.setAtomCache(cache); 
		
		Structure structure = StructureIO.getStructure("1smt");
		
		System.out.println("Testing references in PDB loading with NO alignSeqRes");
		doFullTest(structure);
		
	}
	
	@Test
	public void testCrossReferencesPdbAlignSeqRes() throws IOException, StructureException {
		
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(false);
		
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true); 
		cache.setFileParsingParams(params);
		
		StructureIO.setAtomCache(cache); 
		
		System.out.println("Testing references in PDB loading with alignSeqRes");
		Structure structure = StructureIO.getStructure("1smt");
		
		doFullTest(structure);
		
	}
	
	//@Test
	public void testCrossReferencesRawFile() throws IOException, StructureException {
		// TODO implement
	}

	private void doFullTest(Structure structure) throws StructureException {
		
		System.out.println("Testing references in original structure");
		
		testStructureRefs(structure);
		logger.debug("Original structure mem hashCode: {}",System.identityHashCode(structure));
		
		Structure structureCopy = structure.clone();
		
		logger.debug("Clone structure mem hashCode: {}",System.identityHashCode(structureCopy));
		
		assertNotSame(structure, structureCopy);
		
		System.out.println("Testing references in cloned structure");
		
		testStructureRefs(structureCopy);

		logger.debug("Original structure mem hashCode after cloning: {}",System.identityHashCode(structure));
		
		
		System.out.println("Testing references in original structure after having cloned it");
		// we test again the original after cloning it, perhaps some references were mixed while cloning
		// there is a bug in ChainImpl.clone() that mixes them up!
		testStructureRefs(structure);

		
		System.out.println("Testing references of chain clones");
		for (Chain c:structure.getChains()) {
			Chain clonedChain = (Chain) c.clone();
			testChainRefs(clonedChain);
		}
		
		System.out.println("Testing references in atom arrays");
		for (Chain c:structure.getChains()) {
			Atom[] atomArray = StructureTools.getAllAtomArray(c);
			testAtomArrayRefs(atomArray, c);
		}
		
		
		CrystalBuilder cb = new CrystalBuilder(structure);
		
		StructureInterfaceList interfaces = cb.getUniqueInterfaces();
		
		for (StructureInterface interf:interfaces) {
			testInterfaceRefs(structure, interf);
		}
		
		System.out.println("Testing references in original structure after getUniqueInterfaces");
		testStructureRefs(structure);
	}
	
	private void testStructureRefs(Structure s) throws StructureException {
		
		// structure, chain, group, atom linking
		for (Chain c:s.getChains()) {
			logger.debug("Expected parent structure mem hashCode: {}",System.identityHashCode(s));
			logger.debug("Actual parent structure mem hashCode: {}",System.identityHashCode(c.getParent()));
			assertSame(s, c.getParent()); 
			testChainRefs(c);			
		}
		
		// compounds linking
		for (Compound compound:s.getCompounds()) {
			for (Chain c:compound.getChains()) {
				assertSame(compound, c.getCompound());
				Chain cFromStruc = s.getChainByPDB(c.getChainID());
				assertSame(cFromStruc,c);
			}
		}
	}
	
	private void testChainRefs(Chain c) {
		
		assertNotNull(c.getCompound());
		
		for (Group g:c.getAtomGroups()) {
			assertSame("Failed parent chain for group "+g.toString(), c, g.getChain());
			
			for (Atom a:g.getAtoms()) {
				assertSame(g, a.getGroup());
				assertSame(c, a.getGroup().getChain());
			}
						
			// the SEQRES groups should contain a reference to each and every ATOM group
			// (of course they will also contain some extra groups: those that are only in SEQRES)
			if (c.getSeqResGroups().size()>0) {
				assertTrue("SeqResGroups should contain ATOM group "+g.toString(), 
						containsReference(g, c.getSeqResGroups()) );
			}
		}
		
		for (Group g:c.getSeqResGroups()) {
			assertSame("Failed parent chain for group "+g.toString(), c, g.getChain());
			
			for (Atom a:g.getAtoms()) {
				assertSame(g, a.getGroup());
				assertSame(c, a.getGroup().getChain());
			}
		}
	}

	private void testAtomArrayRefs(Atom[] atoms, Chain c) {
		for (Atom a:atoms) {
			//"Not same: "+c.getChainID()+" ("+c.getAtomLength()+" atomLength, "+c.getSeqResLength()+" seqresLength)"+" | "+a.getGroup().getChain().getChainID()+" ("+a.getGroup().getChain().getAtomLength()+" atomLength, "+a.getGroup().getChain().getSeqResLength()+" seqresLength)"
			assertSame(c, a.getGroup().getChain());
			
		}
		
	}
	
	private void testInterfaceRefs(Structure s, StructureInterface i) throws StructureException {
		
		for (Atom a:i.getMolecules().getFirst()) {
			assertNotNull(a.getGroup().getChain().getCompound());
		}
		
		for (Atom a:i.getMolecules().getSecond()) {
			assertNotNull(a.getGroup().getChain().getCompound());
		}

	}
	
	private static boolean containsReference(Group g, List<Group> list) {
		for (Group group:list) {
			if (group==g) return true;
		}
		return false;
	}
	
}
