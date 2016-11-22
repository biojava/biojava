package org.biojava.nbio.structure.io.mmcif;

import org.junit.Test;
import static org.junit.Assert.*;

import java.io.InputStream;
import java.util.zip.GZIPInputStream;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.io.mmcif.MMcifParser;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifParser;
import org.biojava.nbio.structure.io.mmcif.model.AtomSites;

/**
 * A test for the parsing of _atom_sites category in mmcif, equivalent to SCALE records
 * in PDB files. Currently the category is parsed but is not used to populate any field
 * in the Structure hierarchy. Getting the data from the parser can be useful for things like
 * checking the consistency of the SCALE matrix, see https://github.com/eppic-team/owl/issues/4
 * 
 * 
 * 
 * @author Jose Duarte
 *
 */
public class TestAtomSitesParsing {
	
	@Test
	public void test4hhb() throws Exception {
		
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/4hhb.cif.gz"));
		assertNotNull(inStream);

		MMcifParser mmcifpars = new SimpleMMcifParser();
		SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
		
		mmcifpars.addMMcifConsumer(consumer);

		mmcifpars.parse(inStream) ;
		
		AtomSites atomSites = consumer.getAtomSites();
		
		//System.out.println(consumer.getAtomSites());
		
		Matrix4d m = new Matrix4d(
				Double.parseDouble(atomSites.getFract_transf_matrix11()), Double.parseDouble(atomSites.getFract_transf_matrix12()), Double.parseDouble(atomSites.getFract_transf_matrix13()), Double.parseDouble(atomSites.getFract_transf_vector1()),
				Double.parseDouble(atomSites.getFract_transf_matrix21()), Double.parseDouble(atomSites.getFract_transf_matrix22()), Double.parseDouble(atomSites.getFract_transf_matrix23()), Double.parseDouble(atomSites.getFract_transf_vector2()),
				Double.parseDouble(atomSites.getFract_transf_matrix31()), Double.parseDouble(atomSites.getFract_transf_matrix32()), Double.parseDouble(atomSites.getFract_transf_matrix33()), Double.parseDouble(atomSites.getFract_transf_vector3()),
				0,0,0,1);
		
		// translation components are 0
		assertEquals(0.26656, m.m03, 0.00001);
		assertEquals(0.16413, m.m13, 0.00001);
		assertEquals(0.75059, m.m23, 0.00001);
		
		assertEquals(0.015462, m.m00, 0.00001);
		assertEquals(0.002192, m.m01, 0.00001);
	}

}
