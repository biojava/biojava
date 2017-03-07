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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.contact.StructureInterface;
import org.biojava.nbio.structure.contact.StructureInterfaceList;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.xtal.CrystalBuilder;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class TestStructureCrossReferences {

	private static final Logger logger = LoggerFactory.getLogger(TestStructureCrossReferences.class);

	private static final String PDBCODE1 = "1smt";
	private static final String PDBCODE2 = "2mre";

	@Test
	public void testCrossReferencesMmCif() throws IOException, StructureException {
		boolean emptySeqRes = true;

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);

		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(false); // Store empty seqres groups.
		cache.setFileParsingParams(params);

		StructureIO.setAtomCache(cache);

		Structure structure = StructureIO.getStructure(PDBCODE1);

		//System.out.println("Testing references in mmCIF loading with NO alignSeqRes");
		doFullTest(structure, emptySeqRes);

		structure = StructureIO.getStructure(PDBCODE2); // an NMR entry with 2 chains
		doFullTest(structure, emptySeqRes);

	}

	@Test
	public void testCrossReferencesMmCifAlignSeqRes() throws IOException, StructureException {
		boolean emptySeqRes = false;

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);

		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);

		StructureIO.setAtomCache(cache);

		Structure structure = StructureIO.getStructure(PDBCODE1);

		//System.out.println("Testing references in mmCIF loading with alignSeqRes");
		doFullTest(structure, emptySeqRes);

		structure = StructureIO.getStructure(PDBCODE2); // an NMR entry with 2 chains
		doFullTest(structure, emptySeqRes);


	}


	@Test
	public void testCrossReferencesPdb() throws IOException, StructureException {
		boolean emptySeqRes = true;
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(false);

		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(false);  // Store empty seqres groups
		cache.setFileParsingParams(params);

		StructureIO.setAtomCache(cache);

		Structure structure = StructureIO.getStructure(PDBCODE1);


		doFullTest(structure, emptySeqRes);

		structure = StructureIO.getStructure(PDBCODE2); // an NMR entry with 2 chains
		doFullTest(structure, emptySeqRes);

	}

	@Test
	public void testCrossReferencesPdbAlignSeqRes() throws IOException, StructureException {
		boolean emptySeqRes = false;
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(false);

		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);

		StructureIO.setAtomCache(cache);

		//System.out.println("Testing references in PDB loading with alignSeqRes");
		Structure structure = StructureIO.getStructure(PDBCODE1);

		doFullTest(structure, emptySeqRes);

		structure = StructureIO.getStructure(PDBCODE2); // an NMR entry with 2 chains
		doFullTest(structure, emptySeqRes);


	}

	//@Test
	public void testCrossReferencesRawFile() throws IOException, StructureException {
		// TODO implement
	}

	private void doFullTest(Structure structure, boolean emptySeqRes) throws StructureException {

		//System.out.println("Testing references in original structure");

		testStructureRefs(structure, emptySeqRes);
		logger.debug("Original structure mem hashCode: {}",System.identityHashCode(structure));

		Structure structureCopy = structure.clone();

		logger.debug("Clone structure mem hashCode: {}",System.identityHashCode(structureCopy));

		assertNotSame(structure, structureCopy);

		//System.out.println("Testing references in cloned structure");

		testStructureRefs(structureCopy, emptySeqRes);

		logger.debug("Original structure mem hashCode after cloning: {}",System.identityHashCode(structure));


		//System.out.println("Testing references in original structure after having cloned it");
		// we test again the original after cloning it, perhaps some references were mixed while cloning
		// there is a bug in ChainImpl.clone() that mixes them up!
		testStructureRefs(structure, emptySeqRes);


		//System.out.println("Testing references of chain clones");
		for (Chain c:structure.getChains()) {
			Chain clonedChain = (Chain) c.clone();
			testChainRefs(clonedChain, emptySeqRes);
		}

		//System.out.println("Testing references in atom arrays");
		for (Chain c:structure.getChains()) {
			Atom[] atomArray = StructureTools.getAllAtomArray(c);
			testAtomArrayRefs(atomArray, c);
		}


		CrystalBuilder cb = new CrystalBuilder(structure);

		StructureInterfaceList interfaces = cb.getUniqueInterfaces();

		for (StructureInterface interf:interfaces) {
			testInterfaceRefs(structure, interf);
		}

		//System.out.println("Testing references in original structure after getUniqueInterfaces");
		testStructureRefs(structure, emptySeqRes);
	}

	private void testStructureRefs(Structure s, boolean emptySeqRes) throws StructureException {

		// structure, chain, group, atom linking
		for (Chain c:s.getChains()) {
			logger.debug("Expected parent structure mem hashCode: {}",System.identityHashCode(s));
			logger.debug("Actual parent structure mem hashCode: {}",System.identityHashCode(c.getStructure()));
			assertSame(s, c.getStructure());
			testChainRefs(c, emptySeqRes);
		}

		// compounds linking
		for (EntityInfo compound:s.getEntityInfos()) {
			for (Chain c:compound.getChains()) {
				assertSame(compound, c.getEntityInfo());
				int count = 0;
				for (int modelNr=0;modelNr<s.nrModels();modelNr++) {
					// the chain must be matched by 1 and only 1 chain from all models
					Chain cFromStruc = s.getChain(c.getId(), modelNr);
					if (cFromStruc==c) count++;
				}
				assertEquals("Only 1 chain must match the compound chain for all models",1,count);
			}
		}
	}

	private void testChainRefs(Chain c, boolean emptySeqRes) {

		assertNotNull(c.getEntityInfo());

		for (Group g:c.getAtomGroups()) {
			assertSame("Failed parent chain for group "+g.toString(), c, g.getChain());

			for (Atom a:g.getAtoms()) {
				assertSame(g, a.getGroup());
				assertSame(c, a.getGroup().getChain());
			}

			// the SEQRES groups should contain a reference to each and every ATOM group
			// (of course they will also contain some extra groups: those that are only in SEQRES)
			if (c.getSeqResGroups().size()>0) {
				// we don't want to test hetatoms that are most likely outside of the seqres defined chains
				if (g.getType() == GroupType.HETATM) continue;

				// Only test for expected Atom(s) if we parsed Atoms into Structure.
				if (!emptySeqRes)
				assertTrue("SeqResGroups should contain ATOM group "+g.toString(),
						containsReference(g, c.getSeqResGroups()) );
			}
		}

		for (Group g:c.getSeqResGroups()) {
			// When not aligned, this can be different.
			if (!emptySeqRes)
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
			assertNotNull(a.getGroup().getChain().getEntityInfo());
		}

		for (Atom a:i.getMolecules().getSecond()) {
			assertNotNull(a.getGroup().getChain().getEntityInfo());
		}

	}

	private static boolean containsReference(Group g, List<Group> list) {
		for (Group group:list) {
			if (group==g) return true;
		}
		return false;
	}

}
