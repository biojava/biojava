/*
 * BioJava development code
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
 */
package org.biojava.nbio.structure.cath;

import org.junit.jupiter.api.Test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

public class CathInstallationTest {

	@Test
	public void testParseCathDomainListSuccess() throws IOException {
		String data = "# CATH domain list\n"
				+ "\n"
				+ "1oaiA00     1    10   490    10     1     1     1     1     1   124 1.80\n"
				+ "1oaiA01     1    10   490    10     1     1     1     1     2   150 1.80\n";

		CathInstallation installation = new CathInstallation("");
		BufferedReader reader = new BufferedReader(new StringReader(data));
		installation.parseCathDomainList(reader);

		CathDomain domain = installation.getDomainByCathId("1oaiA00");
		assertNotNull(domain);
		assertEquals("1oaiA00", domain.getDomainName());
		assertEquals(1, domain.getClassId());
		assertEquals(10, domain.getArchitectureId());
		assertEquals(490, domain.getTopologyId());
		assertEquals(10, domain.getHomologyId());
		assertEquals(124, domain.getLength());
		assertEquals(1.80, domain.getResolution(), 0.001);
	}

	@Test
	public void testParseCathDomainListEmptyThrowsException() {
		CathInstallation installation = new CathInstallation("");
		BufferedReader reader = new BufferedReader(new StringReader(""));
		assertThrows(IOException.class, () -> installation.parseCathDomainList(reader));
	}

	@Test
	public void testParseCathDomainListOnlyCommentsAndWhitespaceThrowsException() {
		String data = "# comment 1\n"
				+ "# comment 2\n"
				+ "   \n"
				+ "\t\n";
		CathInstallation installation = new CathInstallation("");
		BufferedReader reader = new BufferedReader(new StringReader(data));
		assertThrows(IOException.class, () -> installation.parseCathDomainList(reader));
	}

	@Test
	public void testParseCathDomainListNoParsableLinesThrowsException() {
		String data = "# comment\n"
				+ "invalid line with too few tokens\n"
				+ "another bad line\n";
		CathInstallation installation = new CathInstallation("");
		BufferedReader reader = new BufferedReader(new StringReader(data));
		IOException exception = assertThrows(IOException.class, () -> installation.parseCathDomainList(reader));
		assertNotNull(exception.getMessage());
	}
}
