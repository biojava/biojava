/**
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

import org.biojava.nbio.structure.scop.RemoteScopInstallation;
import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.scop.ScopInstallation;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Tests {@link ScopInstallation}.
 * @author Spencer Bliven
 * @since 3.0.6
 */
@RunWith(Parameterized.class)
public class RemoteScopInstallationTest extends ScopDatabaseTest {

	public RemoteScopInstallationTest(String tag,ScopDatabase scop) {
		super(tag,scop);
	}

	//@Parameters
	@Parameters(name="{0}")
	public static Collection<Object[]> availableDatabases() {
		ArrayList<Object[]> databases = new ArrayList<Object[]>();
		RemoteScopInstallation scop;
		for(String version : new String[] {
				ScopFactory.LATEST_VERSION,
				ScopFactory.VERSION_1_75A,
				ScopFactory.VERSION_1_75B,
				ScopFactory.VERSION_1_75,
				ScopFactory.VERSION_1_73,
		}) {
			scop = new RemoteScopInstallation();
			scop.setScopVersion(version);
			databases.add(new Object[] {scop.getScopVersion().trim(), scop});
		}
		return databases;
	}
}
