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
 * Created on May 21, 2010
 * Author: Andreas Prlic
 *
 */

package org.biojava.nbio.structure.test.io;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.AfpChainWriter;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.xml.AFPChainXMLParser;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import static org.junit.Assert.assertNotNull;


public class LoadOldXMLfileTest {

	@Test
	public void testLoadOldXMLFile1() throws IOException, StructureException {


		String name1="1P80.D";
		String name2="2IUF.E";

		loadOldXMLFile(name1, name2);


	}

	@Test
	public void testLoadOldXMLFile2() throws IOException, StructureException {


		String name1="1FEZ.A";
		String name2="1O08.A";

		loadOldXMLFile(name1, name2);


	}


	private void loadOldXMLFile(String name1, String name2) throws IOException, StructureException {

		System.out.println("loading " + name1 + " " + name2);
		InputStream inStream = this.getClass().getResourceAsStream("/align/"+name1+"_"+name2+".xml");
		assertNotNull(inStream);

		String xml = convertStreamToString(inStream);

		AtomCache cache = new AtomCache();

		Atom[] ca1 = cache.getAtoms(name1);
		Atom[] ca2 = cache.getAtoms(name2);

		AFPChain afpChain = AFPChainXMLParser.fromXML(xml, ca1, ca2);



		String txt = AfpChainWriter.toWebSiteDisplay(afpChain, ca1, ca2);
		assertNotNull(txt);


	}

	public static String convertStreamToString(InputStream stream) throws IOException {
		BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
		StringBuilder sb = new StringBuilder();

		String line = null;

		while ((line = reader.readLine()) != null) {
		  sb.append(line).append("\n");
		}

		stream.close();


		return sb.toString();
	}
}
