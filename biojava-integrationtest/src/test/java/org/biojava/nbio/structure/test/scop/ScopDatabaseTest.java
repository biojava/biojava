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
 * Created on Sep 30, 2013
 * Author: blivens
 *
 */

package org.biojava.nbio.structure.test.scop;

import org.biojava.nbio.structure.scop.*;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.List;
import java.util.regex.Pattern;

import static org.junit.Assert.*;

/**
 * Generic tests for ScopDatabases. All implementing classes should pass these tests.
 *
 * This abstract class defines the tests for the interface. Each implementation
 * should extend this class and implement the requirements for a parametic test
 * (annotate the class with @RunWith, add an @Parameters static method,
 * and provide a single constructor calling super)
 * @author blivens
 *
 */
@RunWith(Parameterized.class)
public abstract class ScopDatabaseTest {
	protected ScopDatabase scop;
	String tag;

	/**
	 *
	 * @param tag A short string, displayed for failed asserts
	 * @param scop the database instance to test
	 */
	public ScopDatabaseTest(String tag, ScopDatabase scop) {
		if( tag != null) {
			this.tag = "["+tag+"] "; // for messages
		} else {
			this.tag = "";
		}

		this.scop = scop;
	}

	/**
	 * Traverse through the SCOP hierarchy
	 *
	 */
	@Test
	public void traverseHierarchy()
	{
		String pdbId = "4HHB";

		List<ScopDomain> domains = scop.getDomainsForPDB(pdbId);
		assertEquals(tag+"Wrong number of domains",4,domains.size());

		// Check domains (order doesn't matter)
		assertEquals(tag+"Wrong domain","d4hhba_",domains.get(0).getScopId());
		assertEquals(tag+"Wrong domain","d4hhbb_",domains.get(2).getScopId());
		assertEquals(tag+"Wrong domain","d4hhbc_",domains.get(1).getScopId());
		assertEquals(tag+"Wrong domain","d4hhbd_",domains.get(3).getScopId());

		// Check the heirarchy
		ScopNode node = scop.getScopNode(domains.get(0).getSunid());
		ScopDescription desc = scop.getScopDescriptionBySunid(node.getSunid());
		assertEquals(tag,15251,node.getSunid());
		assertEquals(tag,"d4hhba_",desc.getName());
		assertEquals(tag,"4hhb A:",desc.getDescription());
		assertEquals(tag,"a.1.1.2",desc.getClassificationId());

		node = scop.getScopNode(node.getParentSunid());
		desc = scop.getScopDescriptionBySunid(node.getSunid());
		assertEquals(tag,46487,node.getSunid());
		assertEquals(tag,"-",desc.getName());
		assertTrue(tag,Pattern.matches("Human \\(Homo sapiens\\)( \\[TaxId: 9606\\])?",desc.getDescription()));
		assertEquals(tag,"a.1.1.2",desc.getClassificationId());

		node = scop.getScopNode(node.getParentSunid());
		desc = scop.getScopDescriptionBySunid(node.getSunid());
		assertEquals(tag,46486,node.getSunid());
		assertEquals(tag,"-",desc.getName());
		assertEquals(tag,"Hemoglobin, alpha-chain",desc.getDescription());
		assertEquals(tag,"a.1.1.2",desc.getClassificationId());

		node = scop.getScopNode(node.getParentSunid());
		desc = scop.getScopDescriptionBySunid(node.getSunid());
		assertEquals(tag,46463,node.getSunid());
		assertEquals(tag,"-",desc.getName());
		assertEquals(tag,"Globins",desc.getDescription());
		assertEquals(tag,"a.1.1.2",desc.getClassificationId());

		node = scop.getScopNode(node.getParentSunid());
		desc = scop.getScopDescriptionBySunid(node.getSunid());
		assertEquals(tag,46458,node.getSunid());
		assertEquals(tag,"-",desc.getName());
		assertEquals(tag,"Globin-like",desc.getDescription());
		assertEquals(tag,"a.1.1",desc.getClassificationId());


		node = scop.getScopNode(node.getParentSunid());
		desc = scop.getScopDescriptionBySunid(node.getSunid());
		assertEquals(tag,46457,node.getSunid());
		assertEquals(tag,"-",desc.getName());
		assertEquals(tag,"Globin-like",desc.getDescription());
		assertEquals(tag,"a.1",desc.getClassificationId());

		node = scop.getScopNode(node.getParentSunid());
		desc = scop.getScopDescriptionBySunid(node.getSunid());
		assertEquals(tag,46456,node.getSunid());
		assertEquals(tag,"-",desc.getName());
		assertEquals(tag,"All alpha proteins",desc.getDescription());
		assertEquals(tag,"a",desc.getClassificationId());

		// root node
		node = scop.getScopNode(node.getParentSunid());
		desc = scop.getScopDescriptionBySunid(node.getSunid());
		assertEquals(tag,0,node.getSunid());
		assertNull(tag+"Root should not have a description", desc);
	}

	@Test
	public void testNodes() {

		ScopNode node = scop.getScopNode(15251); //4hhb
		assertEquals(tag,15251,node.getSunid());
		assertEquals(tag,46487,node.getParentSunid());

		node = scop.getScopNode(46456); // all-alpha
		assertEquals(tag,46456,node.getSunid());
		assertEquals(tag,0,node.getParentSunid());

		node = scop.getScopNode(0); // root
		assertEquals(tag,0,node.getSunid());
		assertEquals(tag,-1,node.getParentSunid());
		assertEquals(tag+"Wrong number of children",
				// Class I (Artifacts) added in SCOP 2.06
				scop.getScopVersion().compareToIgnoreCase( ScopFactory.VERSION_2_0_6) >= 0 ? 12 : 11,
						node.getChildren().size());

		node = scop.getScopNode(-1); // illegal
		assertNull(tag,node);

		node = scop.getScopNode(Integer.MAX_VALUE); // unused
		assertNull(tag,node);

	}

	/** Get various categories
	 *
	 */
	@Test
	public void testCategories(){
		List<ScopDescription> superfams = scop.getByCategory(ScopCategory.Superfamily);

		assertNotNull(tag,superfams);
		if(scop.getScopVersion().compareToIgnoreCase( ScopFactory.VERSION_1_75) == 0 ) {
			assertEquals(tag,2223,superfams.size());
		} else {
			// defaults for other versions
			assertFalse(tag,superfams.isEmpty());
		}

		List<ScopDescription> folds = scop.getByCategory(ScopCategory.Fold);

		assertNotNull(tag,folds);
		assertFalse(tag,folds.isEmpty());
		if(scop.getScopVersion().compareToIgnoreCase( ScopFactory.VERSION_1_75) == 0 ) {
			assertEquals(tag,1393,folds.size());
		}
	}

	@Test
	public void testComments() {
		// Note that comments change often between SCOP versions.
		// This test is likely to need updating after each SCOPe release.
		List<String> comments;

		// root node
		comments = scop.getComments(0);
		assertTrue(comments.isEmpty());

		//TODO add additional version checks, since comments change a lot

		if(scop.getScopVersion().compareToIgnoreCase( ScopFactory.VERSION_1_71) >= 0 ) {
			comments = scop.getComments(15016);
			assertEquals(tag+"Wrong number of comments", 1, comments.size());
			assertEquals(tag+"Wrong comment", "complexed with cmo, hem", comments.get(0).trim());
			comments = scop.getComments(82738);
			assertEquals(tag+"Wrong number of comments", 1, comments.size());
			assertEquals(tag+"Wrong comment", "inverting reaction mechanism", comments.get(0).trim());
		}
		if(scop.getScopVersion().compareToIgnoreCase( ScopFactory.VERSION_1_73) >= 0 &&
				scop.getScopVersion().compareToIgnoreCase( ScopFactory.VERSION_1_75) <= 0) {
			// Note: only tested so far with 1.75, so may need some modification for earlier versions

			comments = scop.getComments(127355);
			assertEquals(tag+"Wrong number of comments", 2, comments.size());
			assertEquals(tag+"Wrong comment", "automatically matched to d2hbia_", comments.get(0).trim());
			assertTrue(tag+"Wrong comment", Pattern.matches("complexed with hem(; mutant)?", comments.get(1).trim()));
		}
		if(scop.getScopVersion().compareToIgnoreCase( ScopFactory.VERSION_1_75) == 0 ) {
			comments = scop.getComments(160555);
			assertEquals(tag+"Wrong number of comments", 1, comments.size());
			assertEquals(tag+"Wrong comment", "<a href=\"http://pfam.sanger.ac.uk/family?acc=PF06262\">PF06262</a>; DUF1025; minimal zincin fold that retains 3-stranded mixed beta-sheet and contains HExxH motif in the C-terminal helix; no metal ion bound to this motif is observed in the first determined structures", comments.get(0));


		}
		if(scop.getScopVersion().compareToIgnoreCase( ScopFactory.VERSION_1_75B) >= 0 ) {
			// The following were added or modified in 1.75B

			comments = scop.getComments(127355);
			assertEquals(tag+"Wrong number of comments", 2, comments.size());
			assertEquals(tag+"Wrong comment", "automated match to d2hbia_", comments.get(0).trim());
			assertEquals(tag+"Wrong comment", "complexed with hem", comments.get(1).trim());
			// d3ueea_ was added in 1.75B update
			// domain
			comments = scop.getComments(190700);
			assertEquals(tag+"Wrong number of comments", 1, comments.size());
			assertEquals(tag+"Wrong comment", "not a true protein", comments.get(0));

			// fold
			comments = scop.getComments(57923);
			assertEquals(tag+"Wrong number of comments", 1, comments.size());
			assertEquals(tag+"Wrong comment", "metal(zinc)-bound alpha+beta fold", comments.get(0));
		}
	}

}
