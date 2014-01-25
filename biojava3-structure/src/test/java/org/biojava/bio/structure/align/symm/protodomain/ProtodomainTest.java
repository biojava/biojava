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
 * Created on 2012-11-20
 *
 */

package org.biojava.bio.structure.align.symm.protodomain;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureIdentifier;
import org.biojava.bio.structure.align.symm.protodomain.ResourceList.NameProvider;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopFactory;
import org.junit.Before;
import org.junit.Test;

/**
 * A unit test for {@link Protodomain}.
 * @author dmyersturnbull
 *
 */
public class ProtodomainTest {

	private AtomCache cache;

	@Before
	public void setUp() throws Exception {
		ResourceList.set(NameProvider.defaultNameProvider(), ResourceList.DEFAULT_PDB_DIR);
		cache = ResourceList.get().getCache();
	}

	@Test
	public void testConcat() throws Exception {
		AFPChainAndAtoms storedAcaa = ResourceList.get().loadSymm("1ngk.B");
		Protodomain protodomain = Protodomain.fromSymmetryAlignment(storedAcaa.getIdentifier1(), storedAcaa.getAfpChain(), storedAcaa.getCa1(), 2, cache);
	}
	
	@Test
	public void testWholeChain() throws Exception {
		checkAlignSymm("1hiv.A_1-44,A_58-88,A_93-95", "1hiv.A");
	}

	@Test
	public void testWholePdbFile() throws Exception {
		checkAlignSymm("1hiv.A_1-99,B_1-99", "1hiv");
	}

	@Test
	public void testNegativeStart() throws Exception {
		checkAlignSymm("2ags.A_5-20,A_27-57,A_62-87,A_94-114,A_129-139,A_151-165,A_172-294,A_299-373,A_379-391", "d2agsa2"); // issue #3: has insertion codes in the SCOP range; USED TO BE 1qdm.A_3S-37S,A_65S-99S	
	}

	@Test
	public void testInsertionCodesInScopDomain() throws Exception {
		checkAlignSymm("1qdm.A_3S-99S", "d1qdma1"); // issue #3: has insertion codes in the SCOP range; USED TO BE 1qdm.A_3S-37S,A_65S-99S	
	}

	@Test
	public void testHeteroAtomsAtEnd() throws Exception {
		checkAlignSymm("2j8g.A_5-94,A_102-132,A_137-140,A_145-156,A_164-182", "d2j8ga2"); // issue #1: the domain has water atoms at the end of the chain; USED TO BE 2j8g.A_5-170,A_171-182
	}
	
	@Test
	public void testDomainInsertion() throws Exception {
		checkAlignSymm("1w0p.A_221-248,A_257-317,A_328-346,A_544-556,A_564-603,A_615-638,A_648-661,A_669-670,A_676-775", "d1w0pa3"); // issue #2: there is a domain insertion between the start and end
	}

	@Test
	public void testMultiChainDomain() throws Exception {
		checkAlignSymm("1cph.B_4-6,B_11-30,A_1-9", "d1cph.1");
	}

	@Test
	public void testSplicing() throws Exception {
		// with consecutive = 4, this is "3iek.A_17-28,A_56-82,A_87-128,A_133-134,A_150-166,A_173-239,A_244-265,A_273-294,A_320-352,A_369-377"
		AFPChainAndAtoms storedAcaa = ResourceList.get().loadSymm("d3ieka_");
		Protodomain protodomain = Protodomain.fromSymmetryAlignment(storedAcaa.getIdentifier1(), storedAcaa.getAfpChain(), storedAcaa.getCa1(), 1, cache);
		Protodomain spliced5 = protodomain.spliceApproxConsecutive(5);
		assertEquals("3iek.A_17-28,A_56-134,A_150-166,A_173-265,A_273-294,A_320-352,A_369-377", spliced5.toString());
		Protodomain spliced10 = protodomain.spliceApproxConsecutive(10);
		assertEquals("3iek.A_17-28,A_56-134,A_150-294,A_320-352,A_369-377", spliced10.toString());
		Protodomain spliced20 = protodomain.spliceApproxConsecutive(20);
		assertEquals("3iek.A_17-28,A_56-294,A_320-377", spliced20.toString());
		Protodomain spliced100 = protodomain.spliceApproxConsecutive(100);
		assertEquals("3iek.A_17-377", spliced100.toString());
	}
	
	@Test
	public void testSimilarityAlignment() throws Exception {
		checkAlignSim("2bib.A_56-103", "1qdm.A_3S-37S,A_65S-99S", "d2biba2");
	}

	@Test
	public void testSubstructure() throws Exception {
		AFPChainAndAtoms storedAcaa = ResourceList.get().loadSim("1qdm.A_3S-37S,A_65S-99S", "d2biba2");
		Protodomain protodomain = Protodomain.fromReferral(storedAcaa.getIdentifier2(), storedAcaa.getAfpChain(), storedAcaa.getCa2(), cache);
		checkSubstruct("2bib.A_56-79", protodomain, 2, 0); // note that the protodomain is rounded down in all of these
		checkSubstruct("2bib.A_56-71", protodomain, 3, 0);
		checkSubstruct("2bib.A_56-67", protodomain, 4, 0);
		checkSubstruct("2bib.A_56-65", protodomain, 5, 0);
		checkSubstruct("2bib.A_56-63", protodomain, 6, 0);
	}

	@Test
	public void testSubstructureWithIndex() throws Exception {
		/*
		 * Should be:
		 * ngk.B_10-33,B_38-58,B_77-126
		 * ngk.B_10-33,B_38-58,B_77-80
		 * 1ngk.B_81-126
		 */
		AFPChainAndAtoms storedAcaa = ResourceList.get().loadSymm("1ngk.B");
		// order here shouldn't affect the ability to create a substructure, so just say 5
		Protodomain protodomain = Protodomain.fromSymmetryAlignment(storedAcaa.getIdentifier1(), storedAcaa.getAfpChain(), storedAcaa.getCa1(), 5, cache);
		checkSubstruct("1ngk.B_10-33,B_38-58,B_77-80", protodomain, 2, 0);
		checkSubstruct("1ngk.B_81-126", protodomain, 2, 1);
		// TODO check other orders
	}
	
	@Test
	public void testLongGaps() throws Exception {
		checkAlignSymm("3iek.A_17-28,A_56-82,A_87-128,A_133-134,A_150-166,A_173-239,A_244-265,A_273-294,A_320-352,A_369-377", "d3ieka_");
		checkAlignSymm("1uuq.A_40-52,A_58-87,A_96-133,A_143-148,A_171-193,A_199-223,A_230-294,A_301-364,A_369-384,A_395-398,A_412-428", "d1uuqa_");
	}
	
	@Test
	public void testAssorted() throws Exception {
		checkAlignSymm("1u7g.A_13-64,A_77-82,A_101-187,A_194-334", "d1u7ga_");
		checkAlignSymm("1w0p.A_221-248,A_257-317,A_328-346,A_544-556,A_564-603,A_615-638,A_648-661,A_669-670,A_676-775", "d1w0pa3");
		checkAlignSymm("1zgk.A_326-381,A_387-609", "d1zgka1");
		checkAlignSymm("2h6f.B_573-593,B_600-635,B_645-813,B_827-877", "d2h6fb1");
		checkAlignSymm("1juh.A_8-87,A_92-145,A_150-154,A_204-349", "d1juha_");
		checkAlignSymm("1erj.A_335-684,A_691-706", "d1erja_");
	}

	@Test(expected=ProtodomainCreationException.class)
	public void test1Block() throws ProtodomainCreationException {
		AFPChainAndAtoms acaa = ResourceList.get().load(ResourceList.get().openFile("1_block_alignment.xml"), "d1qdma1", "d1qdma1");
		StructureIdentifier identifier = ScopFactory.getSCOP().getDomainByScopID("d1qdma1");
		Protodomain.fromSymmetryAlignment(identifier, acaa.getAfpChain(), acaa.getCa1(), 1, cache); // should throw an exception
	}

	Protodomain checkSubstruct(String string, Protodomain parent, int order, int index) throws ProtodomainCreationException, IOException, StructureException {
//		Protodomain sub = parent.createSubstruct(order, index);
//
//		// check that the protodomain's string is correct
//		assertEquals("The protodomain created has the wrong string", string, sub.toString());
//
//		return sub;
		return null;
	}

	Protodomain checkAlignSim(String string, String name1, String name2) throws StructureException, IOException, ProtodomainCreationException {

		// first we build the protodomain from a stored alignment
		AFPChainAndAtoms storedAcaa = ResourceList.get().loadSim(name1, name2);
		Protodomain protodomain = Protodomain.fromReferral(storedAcaa.getIdentifier2(), storedAcaa.getAfpChain(), storedAcaa.getCa2(), cache);

		// first, check that the protodomain's string is correct
		assertEquals("The protodomain created has the wrong string", string, protodomain.toString());

		// now we make a protodomain from the string (rather than the alignment)
		Protodomain protodomainFromString = Protodomain.fromString(string, storedAcaa.getIdentifier2(), cache);

		// check its string and structure
		assertEquals("The protodomain created from the string has the wrong string", string, protodomainFromString.toString());
		// check the lengths are the same
		assertEquals("The protodomain created from the string has a different length from the protodomain created from the AFPChain", protodomain.getLength(), protodomainFromString.getLength());

		return protodomain;
	}


	Protodomain checkAlignSymm(String string, String scopId) throws ProtodomainCreationException, IOException, StructureException {

		// first we build the protodomain from a stored alignment
		AFPChainAndAtoms storedAcaa = ResourceList.get().loadSymm(scopId);
		Protodomain protodomain = Protodomain.fromSymmetryAlignment(storedAcaa.getIdentifier1(), storedAcaa.getAfpChain(), storedAcaa.getCa1(), 1, cache);
				
		System.out.println(protodomain.toString());
		// first, check that the protodomain's string is correct
		assertEquals("The protodomain created has the wrong string", string, protodomain.toString());

		// now we make a protodomain from the string (rather than the alignment)
		Protodomain protodomainFromString = Protodomain.fromString(string, storedAcaa.getIdentifier1(), cache);

		// check its string and structure
		assertEquals("The protodomain created from the string has the wrong string", string, protodomainFromString.toString());

		// check the lengths are the same
		assertEquals("The protodomain created from the string has a different length from the protodomain created from the AFPChain", protodomain.getLength(), protodomainFromString.getLength());

		return protodomain;
	}

}
