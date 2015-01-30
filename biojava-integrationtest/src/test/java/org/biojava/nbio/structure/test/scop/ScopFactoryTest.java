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
package org.biojava.nbio.structure.test.scop;

import org.biojava.nbio.structure.scop.BerkeleyScopInstallation;
import org.biojava.nbio.structure.scop.RemoteScopInstallation;
import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.junit.Before;
import org.junit.Test;

import java.lang.reflect.Field;
import java.util.HashMap;

import static org.junit.Assert.*;

public class ScopFactoryTest {

	@Before
	public void setUp() throws Exception {
		// reset static values
		Field versionedScopDBs = ScopFactory.class.getDeclaredField("versionedScopDBs");
		versionedScopDBs.setAccessible(true);
		versionedScopDBs.set(null, new HashMap<String, ScopDatabase>());
		Field defaultVersion = ScopFactory.class.getDeclaredField("defaultVersion");
		defaultVersion.setAccessible(true);
		defaultVersion.set(null, ScopFactory.LATEST_VERSION);
	}

	@Test
	public void testVersionCaching() {
		ScopDatabase scop1,scop2,scop3;
		scop1 = ScopFactory.getSCOP(ScopFactory.LATEST_VERSION);
		scop2 = ScopFactory.getSCOP(ScopFactory.LATEST_VERSION);
		assertSame(scop1, scop2);

		scop2 = ScopFactory.getSCOP(ScopFactory.VERSION_1_75);
		assertNotSame(scop1,scop2);
		scop3 = ScopFactory.getSCOP(ScopFactory.VERSION_1_75);
		assertSame(scop2,scop3);

		ScopFactory.setScopDatabase(ScopFactory.VERSION_1_75);
		scop3 = ScopFactory.getSCOP();
		assertSame(scop2,scop3);

		ScopFactory.setScopDatabase(ScopFactory.LATEST_VERSION);
		scop3 = ScopFactory.getSCOP();
		assertSame(scop1,scop3);
	}

	@Test
	public void testVersions() {
		ScopDatabase scop;

		scop = ScopFactory.getSCOP();
		assertEquals(ScopFactory.LATEST_VERSION, scop.getScopVersion());

		scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75);
		assertEquals(ScopFactory.VERSION_1_75, scop.getScopVersion());

		ScopFactory.setScopDatabase(ScopFactory.VERSION_1_75, false);
		scop = ScopFactory.getSCOP();
		assertEquals(ScopFactory.VERSION_1_75, scop.getScopVersion());
		assertSame( RemoteScopInstallation.class,scop.getClass());

		ScopFactory.setScopDatabase(ScopFactory.VERSION_1_75, true);
		scop = ScopFactory.getSCOP();
		assertEquals(ScopFactory.VERSION_1_75, scop.getScopVersion());
		assertSame( BerkeleyScopInstallation.class,scop.getClass());

		ScopFactory.setScopDatabase(ScopFactory.LATEST_VERSION, true);
		scop = ScopFactory.getSCOP();
		assertEquals(ScopFactory.LATEST_VERSION, scop.getScopVersion());
		assertSame( BerkeleyScopInstallation.class,scop.getClass());

	}
}
