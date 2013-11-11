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
package org.biojava.bio.structure.scop;

import java.util.ArrayList;
import java.util.Collection;

import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

/**
 * Tests {@link BerkeleyScopInstallation}.
 * @author dmyerstu
 * @since 3.0.6
 */
@RunWith(Parameterized.class)
public class BerkeleyScopInstallationTest extends ScopDatabaseTest {

	public BerkeleyScopInstallationTest(String tag,ScopDatabase scop) {
		super(tag,scop);
	}
	@Parameters(name="{0}")
	public static Collection<Object[]> availableDatabases() {
		ArrayList<Object[]> databases = new ArrayList<Object[]>();
		ScopInstallation scop;

		for(String version : new String[] {
				//ScopFactory.LATEST_VERSION,
				//ScopFactory.VERSION_1_75A,
				//ScopFactory.VERSION_1_75B,
				ScopFactory.VERSION_1_75,
				ScopFactory.VERSION_1_73,
		}) {
			scop = new BerkeleyScopInstallation();
			scop.setScopVersion(version);
			databases.add(new Object[] {version, scop});
		}
		return databases;
	}
}