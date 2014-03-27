package org.biojava.structure.test.scop;

import static org.junit.Assert.*;

import java.lang.reflect.Field;
import java.util.HashMap;

import org.biojava.bio.structure.scop.BerkeleyScopInstallation;
import org.biojava.bio.structure.scop.RemoteScopInstallation;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopFactory;
import org.junit.Before;
import org.junit.Test;

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
