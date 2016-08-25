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
package org.biojava.nbio.structure.contact;


import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.junit.BeforeClass;
import org.junit.Test;
import org.slf4j.LoggerFactory;

import java.io.IOException;


import static org.junit.Assert.*;




public class TestContactCalc {

	private static final org.slf4j.Logger logger = LoggerFactory.getLogger(TestContactCalc.class);


	private static final String[] INTRACHAIN_TESTSET = {
			"1d2sA",
			"1od3A",
			"1oewA",
			"1w0nA",
			"1wb4A",
			"1wvfA",
			"2gpiA"
			};

	@BeforeClass
	public static void setupBeforeClass() {
		AtomCache cache = new AtomCache();
		FileParsingParameters params = new FileParsingParameters();
		// it is important to set setAlignSeqRes(true), or otherwise contacts won't be stored correctly
		// since the contacts are stored via using SEQRES as identifiers
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);
		StructureIO.setAtomCache(cache);

	}

	@Test
	public void testIntraChainContacts() throws StructureException, IOException {

		String[][] cts = 		{null, {StructureTools.CA_ATOM_NAME} , {StructureTools.CB_ATOM_NAME}};
		double[] cutoffs = 	    { 5.0,   8.0 ,  8.0};

		int[] allCMsizes = new int[INTRACHAIN_TESTSET.length];
		int[] cbCMsizes = new int[INTRACHAIN_TESTSET.length];

		int idx = 0;

		for (String pdbId:INTRACHAIN_TESTSET) {
			logger.info(pdbId+"\t");
			String pdbCode = pdbId.substring(0,4);
			String pdbChainCode = pdbId.substring(4,5);

			Structure structure = StructureIO.getStructure(pdbCode);

			Chain chain = structure.getPolyChainByPDB(pdbChainCode);

			if ( chain!=null) {

				for (int i = 0; i < cts.length; i++) {
					logger.info((cts[i] == null ? "ALL" : cts[i][0]) + "\t" + cutoffs[i] + "\t");

					AtomContactSet atomContacts = null;
					if (cts[i] != null && cts[i][0].equals("CA")) {
						atomContacts = StructureTools.getRepresentativeAtomsInContact(chain, cutoffs[i]);
					} else {
						atomContacts = StructureTools.getAtomsInContact(chain, cts[i], cutoffs[i]);
					}
					GroupContactSet contacts = new GroupContactSet(atomContacts);

					if (cts[i] == null)
						allCMsizes[idx] = contacts.size();

					if (cts[i] != null && cts[i][0].equals("CB"))
						cbCMsizes[idx] = contacts.size();

					int n = chain.getAtomLength();

					assertTrue(contacts.size() > n * 1.20);

					assertTrue(contacts.size() < (n * (n - 1)) / 2);

					// for non-ALL the sizes must be smaller than for ALL
					if (cts[i] != null)
						assertTrue("size for non-ALL contact map (" + contacts.size() + ") should be smaller than for ALL contact map (" + allCMsizes[idx] + ")",
								contacts.size() < allCMsizes[idx]);

					// since the CB contact map will have no contacts for GLYs then the maps should be smaller than the CA
					if (cts[i] != null && cts[i][0].equals("CA"))
						assertTrue("size for CA contact map (" + contacts.size() + ") should be larger than for CB contact map (" + cbCMsizes[idx] + ")",
								contacts.size() > cbCMsizes[idx]);
				}
			}
			logger.info("");
			idx++;
		}

	}

	@Test
	public void testInterChainContacts3HBX() throws StructureException, IOException {

		// 3 interfaces in the AU are NCS equivalent, they should have similar numbers of contacts
		Structure structure = StructureIO.getStructure("3hbx");

		AtomContactSet atomContacts1 = StructureTools.getAtomsInContact(structure.getPolyChainByPDB("A"), structure.getPolyChainByPDB("B"), 5.5, false);
		AtomContactSet atomContacts2 = StructureTools.getAtomsInContact(structure.getPolyChainByPDB("E"), structure.getPolyChainByPDB("F"), 5.5, false);
		AtomContactSet atomContacts3 = StructureTools.getAtomsInContact(structure.getPolyChainByPDB("C"), structure.getPolyChainByPDB("D"), 5.5, false);

		logger.info("AU interfaces of 3hbx, number of atom contacts: "+atomContacts1.size()+", "+atomContacts2.size()+", "+atomContacts3.size());

		assertTrue(Math.abs(atomContacts1.size()-atomContacts2.size())<40);
		assertTrue(Math.abs(atomContacts1.size()-atomContacts3.size())<40);
		assertTrue(Math.abs(atomContacts2.size()-atomContacts3.size())<40);

		GroupContactSet contacts1 = new GroupContactSet(atomContacts1);
		GroupContactSet contacts2 = new GroupContactSet(atomContacts2);
		GroupContactSet contacts3 = new GroupContactSet(atomContacts3);

		logger.info("AU interfaces of 3hbx, number of residue contacts: "+contacts1.size()+", "+contacts2.size()+", "+contacts3.size());

		assertTrue(Math.abs(contacts1.size()-contacts2.size())<10);
		assertTrue(Math.abs(contacts1.size()-contacts3.size())<10);
		assertTrue(Math.abs(contacts2.size()-contacts3.size())<10);

		assertTrue(contacts1.size()>1);
		assertTrue(contacts2.size()>1);
		assertTrue(contacts3.size()>1);

		assertTrue(contacts1.size()<atomContacts1.size()/10);
		assertTrue(contacts2.size()<atomContacts2.size()/10);
		assertTrue(contacts3.size()<atomContacts3.size()/10);

	}

	@Test
	public void testIntraChainContactsVsDistMatrix1SMT() throws IOException, StructureException {

		double cutoff = 5;

		Structure structure = StructureIO.getStructure("1smt");

		Chain chain = structure.getPolyChainByPDB("A");

		logger.info("Intra-chain contacts calculation vs distance matrix for 1smtA");

		checkContactsVsDistMatrix(chain, cutoff);
	}

	@Test
	public void testIntraChainContactsVsDistMatrix2TRX() throws IOException, StructureException {

		double cutoff = 5;

		Structure structure = StructureIO.getStructure("2trx");

		Chain chain = structure.getPolyChainByPDB("A");

		logger.info("Intra-chain contacts calculation vs distance matrix for 2trxA");

		checkContactsVsDistMatrix(chain, cutoff);
	}

	@Test
	public void testIntraChainContactsVsDistMatrix1SU4() throws IOException, StructureException {

		double cutoff = 5;

		Structure structure = StructureIO.getStructure("1su4");

		Chain chain = structure.getPolyChainByPDB("A");

		logger.info("Intra-chain contacts calculation vs distance matrix for 1su4A");

		checkContactsVsDistMatrix(chain, cutoff);
	}

	private void checkContactsVsDistMatrix(Chain chain, double cutoff) {
		long start = System.currentTimeMillis();
		AtomContactSet atomContacts = StructureTools.getAtomsInContact(chain, cutoff);
		long end = System.currentTimeMillis();
		System.out.printf("Calculated contacts in %.3f s\n",((end-start)/1000.0));

		start = System.currentTimeMillis();
		Atom[] atoms = StructureTools.getAllNonHAtomArray(chain,false);
		double[][] distMatrix = calcDistanceMatrix(atoms);
		end = System.currentTimeMillis();
		System.out.printf("Calculated distance matrix in %.3f s\n",((end-start)/1000.0));

		logger.info("(number of atoms: "+atoms.length+")");

		for (int i=0;i<atoms.length;i++) {
			for (int j=i+1;j<atoms.length;j++) {
				if (distMatrix[i][j]<cutoff) {

					assertTrue(atomContacts.hasContact(atoms[i],atoms[j]));

					assertEquals(distMatrix[i][j], atomContacts.getContact(atoms[i],atoms[j]).getDistance(), 0.00001);

				} else {
					assertFalse(atomContacts.hasContact(atoms[i],atoms[j]));
				}
			}
		}

	}

	@Test
	public void testInterChainContactsVsDistMatrix2TRX() throws IOException, StructureException {

		double cutoff = 5;

		Structure structure = StructureIO.getStructure("2trx");

		logger.info("Inter-chain contacts calculation vs distance matrix for 2trx A-B");

		Chain chain1 = structure.getPolyChainByPDB("A");
		Chain chain2 = structure.getPolyChainByPDB("B");

		long start = System.currentTimeMillis();
		AtomContactSet atomContacts = StructureTools.getAtomsInContact(chain1, chain2, cutoff, false);
		long end = System.currentTimeMillis();
		System.out.printf("Calculated contacts in %.3f s\n",((end-start)/1000.0));

		start = System.currentTimeMillis();
		Atom[] atoms1 = StructureTools.getAllNonHAtomArray(chain1,false);
		Atom[] atoms2 = StructureTools.getAllNonHAtomArray(chain2,false);
		double[][] distMatrix = calcDistanceMatrix(atoms1, atoms2);
		end = System.currentTimeMillis();
		System.out.printf("Calculated distance matrix in %.3f s\n",((end-start)/1000.0));

		logger.info("(number of atoms: "+atoms1.length+", "+atoms2.length+")");

		for (int i=0;i<atoms1.length;i++) {
			for (int j=0;j<atoms2.length;j++) {
				if (distMatrix[i][j]<cutoff) {

					assertTrue(atomContacts.hasContact(atoms1[i],atoms2[j]));

					assertEquals(distMatrix[i][j], atomContacts.getContact(atoms1[i],atoms2[j]).getDistance(), 0.00001);

				} else {
					assertFalse(atomContacts.hasContact(atoms1[i],atoms2[j]));
				}
			}
		}
	}

	private double[][] calcDistanceMatrix(Atom[] atoms) {

		double[][] distMatrix = new double[atoms.length][atoms.length];

		for (int i=0;i<atoms.length;i++) {
			for (int j=i+1;j<atoms.length;j++) {
				distMatrix[i][j] = Calc.getDistance(atoms[i], atoms[j]);
			}
		}
		return distMatrix;
	}

	private double[][] calcDistanceMatrix(Atom[] atoms1, Atom[] atoms2) {

		double[][] distMatrix = new double[atoms1.length][atoms2.length];

		for (int i=0;i<atoms1.length;i++) {
			for (int j=0;j<atoms2.length;j++) {
				distMatrix[i][j] = Calc.getDistance(atoms1[i], atoms2[j]);
			}
		}
		return distMatrix;
	}
}
