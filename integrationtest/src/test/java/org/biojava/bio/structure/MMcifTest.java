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
 * created at Apr 26, 2008
 */
package org.biojava.bio.structure;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;

import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileParser;
import org.biojava.bio.structure.io.mmcif.MMcifParser;
import org.biojava.bio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.bio.structure.io.mmcif.SimpleMMcifParser;

import junit.framework.TestCase;

public class MMcifTest extends TestCase {

	protected boolean headerOnly;

	public MMcifTest(){
		super();
		setHeaderOnly(false);
	}



	public boolean isHeaderOnly() {
		return headerOnly;
	}



	public void setHeaderOnly(boolean headerOnly) {
		this.headerOnly = headerOnly;
	}

	public void testLoad(){
		
		// a structure with microheterogeneity
		//comparePDB2cif("2CI1","A");
		
		// test a simple protein
		comparePDB2cif("5pti","A");

		// test a protein with modified residues
		comparePDB2cif("1a4w","L");
		comparePDB2cif("1a4w","H");
		comparePDB2cif("1a4w","I");
                //non-standard encoded amino acid
                comparePDB2cif("1fdo","A");

		// test a DNA binding protein
		comparePDB2cif("1j59","A");
		//comparePDB2cif("1j59","B");
		//comparePDB2cif("1j59","C");
		//comparePDB2cif("1j59","D");
		comparePDB2cif("1j59","E");
		//comparePDB2cif("1j59","F");

		// test a NMR protein
		comparePDB2cif("2kc9","A");

		
	}


	protected void comparePDB2cif(String id, String chainId){
		String fileName = "/"+id+".cif";
		InputStream inStream = this.getClass().getResourceAsStream(fileName);
		assertNotNull("Could not find file " + fileName + ". Config problem?" , inStream);

		MMcifParser parser = new SimpleMMcifParser();

		SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
		FileParsingParameters params = new FileParsingParameters();
		params.setHeaderOnly(headerOnly);
		consumer.setFileParsingParameters(params);
		
		
		parser.addMMcifConsumer(consumer);
		try {
			parser.parse(new BufferedReader(new InputStreamReader(inStream)));
		} catch (IOException e){
			fail(e.getMessage());
		}
		// remove to avoid memory leaks
		parser.clearConsumers();
		Structure cifStructure = consumer.getStructure();
		assertNotNull(cifStructure);


		// load the PDB file via the PDB parser
		Structure pdbStructure = null;
		InputStream pinStream = this.getClass().getResourceAsStream("/"+id+".pdb");
		assertNotNull(inStream);

		PDBFileParser pdbpars = new PDBFileParser();
		pdbpars.setFileParsingParameters(params);
		
		try {
			pdbStructure = pdbpars.parsePDBFile(pinStream) ;
		} catch (IOException e) {
			e.printStackTrace();
		}

		assertNotNull(pdbStructure);

		// now compare the results
		try {

			// chech NMR data
			assertEquals(id + ": the isNMR flag is not the same!", pdbStructure.isNmr(), cifStructure.isNmr());

			if ( pdbStructure.isNmr()){
				assertEquals(id + ": the nr of NMR models is not the same!", pdbStructure.nrModels(), pdbStructure.nrModels());
				checkNMR(pdbStructure);
				checkNMR(cifStructure);
			}

			//System.out.println(pdbStructure);
			//System.out.println(cifStructure);

			// compare amino acids in chain 1:
			Chain a_pdb = pdbStructure.getChainByPDB(chainId);
			Chain a_cif = cifStructure.getChainByPDB(chainId);
			//System.out.println(a_pdb.getAtomGroups());
			
				//System.out.println(id + "_" + chainId + " pdb atom groups: " + a_pdb.getAtomGroups(GroupType.AMINOACID).size());
				//System.out.println(id + "_" + chainId + " cif atom groups: " + a_cif.getAtomGroups(GroupType.AMINOACID).size());
			
			//for (Group g: a_cif.getAtomGroups()){
			//	System.out.println(g);
			//}
			//System.out.println("--");
			String pdb_SEQseq = a_pdb.getSeqResSequence();
			
			String cif_SEQseq = a_cif.getSeqResSequence();

//                        System.out.println(id + "_" + chainId + " pdbSEQ: " + pdb_SEQseq);
//                        System.out.println(id + "_" + chainId + " cifSEQ: " + cif_SEQseq);
			
			assertEquals(id + ": the SEQRES sequences don't match!", pdb_SEQseq,cif_SEQseq);
			
			assertEquals(id + ":  The nr of ATOM groups does not match!",a_pdb.getAtomGroups(GroupType.AMINOACID).size(),a_cif.getAtomGroups(GroupType.AMINOACID).size()  );
			
			// actually this check not necessarily works, since there can be waters in PDB that we don;t deal with yet in cif...
			//assertEquals("the nr of ATOM record groups is not the same!" , a_pdb.getAtomLength(),a_cif.getAtomLength());
			for (int i = 0 ; i < a_pdb.getAtomGroups(GroupType.AMINOACID).size(); i++){				
				Group gp = a_pdb.getAtomGroups(GroupType.AMINOACID).get(i);
				
				List<Group> cifGroups = a_cif.getAtomGroups(GroupType.AMINOACID);					
				Group gc = cifGroups.get(i);
				checkGroups(gp, gc);
			}



			String pdb_seq = a_pdb.getAtomSequence();
			String cif_seq = a_cif.getAtomSequence();

			//System.out.println(pdb_seq);
			//System.out.println(cif_seq);

			assertEquals("the sequences obtained from PDB and mmCif don't match!",pdb_seq, cif_seq);

			List<DBRef> pdb_dbrefs= pdbStructure.getDBRefs();
			List<DBRef> cif_dbrefs= cifStructure.getDBRefs();

			assertEquals("nr of DBrefs found does not match!", pdb_dbrefs.size(),cif_dbrefs.size());

			DBRef p = pdb_dbrefs.get(0);
			DBRef c = cif_dbrefs.get(0);

			//System.out.println(p.toPDB());
			//System.out.println(c.toPDB());
			String pdb_dbref = p.toPDB();
			String cif_dbref = c.toPDB();
			assertEquals("DBRef is not equal",pdb_dbref,cif_dbref);

			PDBHeader h1 = pdbStructure.getPDBHeader();
			PDBHeader h2 = cifStructure.getPDBHeader();

			//compareString(h1.toPDB() ,h2.toPDB());
			//System.out.println(h1.toPDB());
			//System.out.println(h2.toPDB());
			if ( ! h1.toPDB().toUpperCase().equals(h2.toPDB().toUpperCase()) ){
				System.err.println(h1.toPDB());
				System.err.println(h2.toPDB());
				compareString(h1.toPDB(), h2.toPDB());
			}
			assertEquals("the PDBHeader.toPDB representation is not equivalent", h1.toPDB().toUpperCase(),h2.toPDB().toUpperCase());

			// and the ultimate test!
			// but we are not there yet...
			// TODO: still need to parse SSBOND equivalent info from cif files...
			//assertEquals("the Structure.toPDB representation is not equivalent", pdbStructure.toPDB(),cifStructure.toPDB());

		} catch (StructureException ex){
			fail(ex.getMessage() + " for PDB: " + id);
		}

	}

	private void checkGroups(Group g1, Group g2){

		//System.out.print("comparing " +g1 + " " + g2);
		String pdbId1 = g1.getChain().getParent().getPDBCode();
		String pdbId2 = g1.getChain().getParent().getPDBCode();
		assertEquals(pdbId1, pdbId2);
		
		assertEquals(g1.getType(),g2.getType());
		assertEquals(g1.getResidueNumber().getSeqNum(),g2.getResidueNumber().getSeqNum());
		assertEquals(g1.getResidueNumber().getInsCode(),g2.getResidueNumber().getInsCode());
		assertEquals(g1.getPDBName(),g2.getPDBName());
		assertEquals(g1.has3D(),g2.has3D());
		
		assertEquals(g1.hasAltLoc(), g2.hasAltLoc());
		
		assertEquals(pdbId1 + ":" + g1 + " - " + pdbId2+":"+ g2,g1.getAltLocs().size(), g2.getAltLocs().size());
		
		assertEquals(pdbId1 + ":" + g1 + " - " + pdbId2+":"+ g2 , g1.getAtoms().size(), g2.getAtoms().size());
		if ( g1.has3D()){
			try {
				Atom a1 = g1.getAtom(0);
				Atom a2 = g2.getAtom(0);
				assertEquals(a1.getX(),a2.getX());
				assertEquals(a1.getOccupancy(),a2.getOccupancy());
				assertEquals(a1.getTempFactor(),a2.getTempFactor());
				assertEquals(a1.getFullName(),a2.getFullName());
			} catch (StructureException e){
				fail(e.getMessage());
			}

		}
		//System.out.println(" ... done");

	}

	private void checkNMR(Structure s){
		assertTrue(s.isNmr());

		int models = s.nrModels();
		assertTrue(models > 0);

		List<Chain> model0 = s.getModel(0);

		// compare with all others
		for (int i = 1 ; i < models; i++){
			List<Chain> modelX = s.getModel(i);

			assertEquals(model0.size(),modelX.size());

			// compare lengths:
			for (int j=0 ; j< model0.size();j++){
				Chain c1 = model0.get(j);
				Chain cx = modelX.get(j);
				assertEquals(c1.getAtomLength(),cx.getAtomLength());
				// can;t compare seq res, since this is only done for 1st...
				//assertEquals("c1.getSeqResLength(),cx.getSeqResLength());
				assertEquals(c1.getAtomSequence(),cx.getAtomSequence());
				assertEquals(c1.getAtomGroups("amino").size(),cx.getAtomGroups("amino").size());
				assertEquals(c1.getAtomGroups(GroupType.AMINOACID).size(),cx.getAtomGroups(GroupType.AMINOACID).size());
				assertEquals(c1.getAtomGroups(GroupType.NUCLEOTIDE).size(),cx.getAtomGroups(GroupType.NUCLEOTIDE).size());
				assertEquals(c1.getAtomGroups(GroupType.HETATM).size(),cx.getAtomGroups(GroupType.HETATM).size());
			}


		}
	}

	private void compareString(String t, String pdb){
		for (int i =0 ; i < t.length() ; i++){
			System.out.println(">"+t.charAt(i)+":"+ pdb.charAt(i)+"<");
			if ( Character.toUpperCase(t.charAt(i)) != Character.toUpperCase(pdb.charAt(i))){
				System.out.println("Mismatch!");
				break;
			}
		}
	}


}
