/*
 *			BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *	  http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *	  http://www.biojava.org/
 *
 * Created on Jun 2, 2010
 * Author: Jianjiong Gao
 *
 */

package org.biojava.nbio.protmod.structure;

import org.biojava.nbio.protmod.*;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 *
 * @author Jianjiong Gao
 * @since 3.0
 */
public class ProteinModificationRegistryTest {

	//private static final Logger logger = LoggerFactory.getLogger(ProteinModificationRegistryTest.class);

	/**
	 * Note: if you change this unit test, also change the cook book:
	 * http://www.biojava.org/wiki/BioJava:CookBook3:AddProtMod
	 */
	@Test
	public void testRegisterModification() {
		// define the involved components, in this case two cystines (CYS)
		List<Component> components = new ArrayList<Component>(2);
		components.add(Component.of("CYS"));
		components.add(Component.of("CYS"));

		// define the atom linkages between the components, in this case the SG atoms on both CYS groups
		ModificationLinkage linkage = new ModificationLinkage(components, 0, "SG", 1, "SG");

		// define the modification condition, i.e. what components are involved and what atoms are linked between them
		ModificationCondition condition = new ModificationConditionImpl(components, Collections.singletonList(linkage));

		// build a modification
		ProteinModification mod =
			new ProteinModificationImpl.Builder("0018_test",
				ModificationCategory.CROSS_LINK_2,
				ModificationOccurrenceType.NATURAL,
				condition)
			.setDescription("A protein modification that effectively cross-links two L-cysteine residues to form L-cystine.")
			.setFormula("C 6 H 8 N 2 O 2 S 2")
			.setResidId("AA0025")
			.setResidName("L-cystine")
			.setPsimodId("MOD:00034")
			.setPsimodName("L-cystine (cross-link)")
			.setSystematicName("(R,R)-3,3'-disulfane-1,2-diylbis(2-aminopropanoic acid)")
			.addKeyword("disulfide bond")
			.addKeyword("redox-active center")
			.build();

		//register the modification
		ProteinModificationRegistry.register(mod);
		Assert.assertNotNull(ProteinModificationRegistry.getById("0018_test"));
	}

	/**
	 * Test the initialization registry of common protein modifications.
	 * Note: if you change this unit test, also change the cook book:
		 * http://www.biojava.org/wiki/BioJava:CookBook3:SupportedProtMod
	 */
	@Test
	public void testRegisterCommonModification() {
		Set<ProteinModification> mods = ProteinModificationRegistry.allModifications();
		Assert.assertTrue(mods != null && !mods.isEmpty());

//		logger.info("There are totally "+mods.size()
//				+" protein modifications registered.");
//
//		printModifications(mods);
	}

	/**
	 * Note: if you change this unit test, also change the cook book:
	 * http://www.biojava.org/wiki/BioJava:CookBook3:SupportedProtMod
	 */
	@Test
	public void testGetBy() {
		ProteinModification mod;
		Set<ProteinModification> mods;

		mod = ProteinModificationRegistry.getById("0001");
		Assert.assertNotNull(mod);

		// a set of protein modifications by RESID ID
		mods = ProteinModificationRegistry.getByResidId("AA0151");
		Assert.assertNotNull(mods);

		// a set of protein modifications by PSI-MOD ID
		mods = ProteinModificationRegistry.getByPsimodId("MOD:00305");
		Assert.assertNotNull(mods);

		// a set of protein modifications by PDBCC ID
		mods = ProteinModificationRegistry.getByPdbccId("SEP");
		Assert.assertNotNull(mods);

		// a set of protein modifications by category
		mods = ProteinModificationRegistry.getByCategory(ModificationCategory.ATTACHMENT);
		Assert.assertNotNull(mods);

		// a set of protein modifications by occurrence type
		mods = ProteinModificationRegistry.getByOccurrenceType(ModificationOccurrenceType.NATURAL);
		Assert.assertNotNull(mods);

		// a set of protein modifications by a keyword
		mods = ProteinModificationRegistry.getByKeyword("phosphoprotein");
		Assert.assertNotNull(mods);

		// a set of protein modifications by involved components.
		mods = ProteinModificationRegistry.getByComponent(Component.of("FAD"));
		Assert.assertNotNull(mods);
	}

//	/**
//	 * Print modifications.
//	 * @param mods {@link ProteinModification}s.
//	 */
//	private void printModifications(Set<ProteinModification> mods) {
//		for (ProteinModification mod:mods) {
//			logger.info(mod);
//		}
//	}

}
