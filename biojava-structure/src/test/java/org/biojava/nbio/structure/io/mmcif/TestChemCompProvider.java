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
package org.biojava.nbio.structure.io.mmcif;

import static org.junit.Assert.assertEquals;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class TestChemCompProvider {
	private static final Logger s_logger = LoggerFactory.getLogger(TestChemCompProvider.class);

	// Short test with bad ligand name (QNA is bogus)
	final	String DNAexample =

	"ATOM      1  H   MET A   1      11.756 -15.759  11.647  1.00  7.95\n" +
	"ATOM      2  N   MET A   1      12.461 -16.373  11.329  1.00  7.95\n" +
	"ATOM      3  CA  MET A   1      12.297 -17.782  11.674  1.00  7.95\n" +
	"ATOM      4  C   MET A   1      10.999 -18.339  11.099  1.00  7.95\n" +
	"ATOM      5  O   MET A   1      10.716 -18.148   9.918  1.00  7.95\n" +
	"ATOM      6  CB  MET A   1      12.314 -17.962  13.194  1.00  7.95\n" +
	"ATOM      7  CG  MET A   1      13.685 -17.628  13.781  1.00  7.95\n" +
	"ATOM      8  SD  MET A   1      13.679 -17.742  15.584  1.00  7.95\n" +
	"ATOM      9  CE  MET A   1      12.642 -16.314  15.948  1.00  7.95\n" +
	"ATOM     10  H   GLN A   2      10.381 -19.508  12.646  1.00  6.82\n" +
	"ATOM     11  N   GLN A   2      10.134 -19.027  11.820  1.00  6.82\n" +
	"ATOM     12  CA  GLN A   2       8.782 -19.063  11.370  1.00  6.82\n" +
	"ATOM     13  C   GLN A   2       7.791 -18.046  11.927  1.00  6.82\n" +
	"ATOM     14  O   GLN A   2       6.616 -18.069  11.565  1.00  6.82\n" +
	"ATOM     15  CB  GLN A   2       8.285 -20.485  11.642  1.00  6.82\n" +
	"ATOM     16  CG  GLN A   2       8.996 -21.508  10.755  1.00  6.82\n" +
	"ATOM     17  CD  GLN A   2       8.708 -21.249   9.280  1.00  6.82\n" +
	"ATOM     18  NE2 GLN A   2       9.730 -21.006   8.488  1.00  6.82\n" +
	"ATOM     19  OE1 GLN A   2       7.563 -21.267   8.850  1.00  6.82\n" +
	"ATOM     20  H   MET A   3       9.163 -17.424  13.270  1.00  5.82\n" +
	"ATOM     21  N   MET A   3       8.334 -17.182  12.793  1.00  5.82\n" +
	"ATOM     22  CA  MET A   3       7.753 -15.950  13.031  1.00  5.82\n" +
	"ATOM     23  C   MET A   3       8.145 -14.729  12.164  1.00  5.82\n" +
	"ATOM     24  O   MET A   3       7.643 -13.630  12.384  1.00  5.82\n" +
	"ATOM     25  CB  MET A   3       8.025 -15.652  14.507  1.00  5.82\n" +
	"ATOM     26  CG  MET A   3       7.282 -16.626  15.421  1.00  5.82\n"+
	"ATOM     27  SD  MET A   3       5.496 -16.545  15.162  1.00  5.82\n"+
	"ATOM     28  CE  MET A   3       5.189 -14.910  15.855  1.00  5.82\n" +
	"TER\n" +
	"HETATM  101  O5' QNA A  1      15.062 -23.351   2.519  1.00 66.98           O\n" +
	"HETATM  102  C5' QNA A  1      14.372 -23.705   1.300  1.00 68.27           C\n" +
	"HETATM  103  C4' QNA A  1      14.836 -22.832   0.142  1.00 68.36           C\n" +
	"HETATM  104  O4' QNA A  1      14.402 -23.357  -1.153  1.00 68.93           O\n" +
	"HETATM  105  C3' QNA A  1      14.235 -21.444   0.256  1.00 68.53           C\n" +
	"HETATM  106  O3' QNA A  1      15.060 -20.416  -0.271  1.00 66.63           O\n" +
	"HETATM  107  C2' QNA A  1      12.938 -21.603  -0.522  1.00 69.24           C\n" +
	"HETATM  108  C1' QNA A  1      13.336 -22.562  -1.653  1.00 69.84           C\n" +
	"HETATM  109  N1  QNA A  1      12.146 -23.404  -2.068  1.00 70.20           N\n" +
	"HETATM  110  C2  QNA A  1      11.960 -23.752  -3.392  1.00 69.91           C\n" +
	"HETATM  111  O2  QNA A  1      12.771 -23.504  -4.264  1.00 69.83           O\n" +
	"HETATM  112  N3  QNA A  1      10.812 -24.471  -3.645  1.00 69.68           N\n" +
	"HETATM  113  C4  QNA A  1       9.831 -24.834  -2.730  1.00 70.01           C\n" +
	"HETATM  114  O4  QNA A  1       8.846 -25.491  -3.072  1.00 70.53           O\n" +
	"HETATM  115  C5  QNA A  1      10.079 -24.409  -1.385  1.00 70.05           C\n" +
	"HETATM  116  C7  QNA A  1       9.085 -24.740  -0.319  1.00 69.81           C\n" +
	"HETATM  117  C6  QNA A  1      11.190 -23.716  -1.124  1.00 70.41           C\n" +
	"HETATM  118  P   QNA A  2      14.663 -18.896   0.143  1.00 67.23           P\n" +
	"HETATM  119  OP1 QNA A  2      15.930 -18.124   0.315  1.00 67.71           O\n" +
	"HETATM  120  OP2 QNA A  2      13.607 -18.874   1.201  1.00 65.55           O\n" +
	"HETATM  121  O5' QNA A  2      13.985 -18.327  -1.199  1.00 62.80           O\n" +
	"HETATM  122  C5' QNA A  2      14.827 -18.399  -2.335  1.00 53.86           C\n" +
	"HETATM  123  C4' QNA A  2      13.982 -18.288  -3.576  1.00 47.13           C\n" +
	"HETATM  124  O4' QNA A  2      13.020 -19.368  -3.719  1.00 44.32           O\n" +
	"HETATM  125  C3' QNA A  2      13.193 -17.012  -3.509  1.00 43.63           C\n" +
	"HETATM  126  O3' QNA A  2      13.373 -16.484  -4.756  1.00 42.30           O\n" +
	"HETATM  127  C2' QNA A  2      11.765 -17.470  -3.280  1.00 41.55           C\n" +
	"HETATM  128  C1' QNA A  2      11.759 -18.797  -4.025  1.00 40.38           C\n" +
	"HETATM  129  N1  QNA A  2      10.633 -19.685  -3.585  1.00 35.19           N\n" +
	"HETATM  130  C2  QNA A  2       9.891 -20.348  -4.563  1.00 31.53           C\n" +
	"HETATM  131  O2  QNA A  2      10.201 -20.178  -5.738  1.00 32.59           O\n" +
	"HETATM  132  N3  QNA A  2       8.861 -21.155  -4.207  1.00 28.77           N\n" +
	"HETATM  133  C4  QNA A  2       8.559 -21.301  -2.925  1.00 29.06           C\n" +
	"HETATM  134  N4  QNA A  2       7.535 -22.097  -2.627  1.00 27.42           N\n" +
	"HETATM  135  C5  QNA A  2       9.289 -20.623  -1.901  1.00 31.45           C\n" +
	"HETATM  136  C6  QNA A  2      10.305 -19.829  -2.271  1.00 33.51           C\n";

	@Test
	public void testZipChemCompProvider() throws IOException {
		// Test file input from a stream created from a string.
		InputStream testPDB = new ByteArrayInputStream(DNAexample.getBytes());

		// we just need this to track where to store PDB files
		// this checks the PDB_DIR property (and uses a tmp location if not set)
		UserConfiguration config = new UserConfiguration();
		String cachePath = config.getCacheFilePath();

		// Setup a ChemCompProvider
		Path pdbdir = Paths.get(cachePath);
		Path chemComp = pdbdir.resolve("chemcomp.zip");

		System.out.println("Using PDB_DIR=" + pdbdir.toString() + " as temporary file directory");
		ZipChemCompProvider zp = new ZipChemCompProvider(chemComp.toString(), pdbdir.toString());
		// Keep the .cif.gz files - avoid re-downloading for unit testing.
		zp.setRemoveCif(false);
		ChemCompGroupFactory.setChemCompProvider(zp);

		// Parameters
		FileParsingParameters params = new FileParsingParameters();

		PDBFileReader pdbreader = new PDBFileReader();
		pdbreader.setFileParsingParameters(params);

		Structure s = pdbreader.getStructure(testPDB);
		assertEquals(3, s.getChainByIndex(0).getAtomGroups().size());
		// Not wanted here for testing, but useful for cleaning up downloaded .cif.gz files.
		// ZipChemCompProvider.purgeTempFiles(pdbdir.toString());
	}

	@Test
	public void testNormalStructure() throws StructureException, IOException {
		// we just need this to track where to store PDB files
		// this checks the PDB_DIR property (and uses a tmp location if not set)
		UserConfiguration config = new UserConfiguration();
		String cachePath = config.getCacheFilePath();

		// Setup a ChemCompProvider
		Path pdbdir = Paths.get(cachePath);
		Path chemComp = pdbdir.resolve("chemcomp.zip");

		System.out.println("Using PDB_DIR=" + pdbdir.toString() + " as temporary file directory");

		AtomCache cache = new AtomCache();
		StructureIO.setAtomCache(cache);
		FileParsingParameters params = cache.getFileParsingParams();
		params.setParseBioAssembly(true);

		/*
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());
		StructureIO.setAtomCache(cache);
		cache.setUseMmCif(true);
		long startTime = System.currentTimeMillis();
		Structure sCif = StructureIO.getStructure("4HHM");
		long finishTime = System.currentTimeMillis();
		s_logger.info("DownloadChemComp time: "+(finishTime-startTime)+ " ms");
		*/

		ZipChemCompProvider zp = new ZipChemCompProvider(chemComp.toString(), pdbdir.toString());
		// Keep the .cif.gz files - avoid re-downloading for unit testing.
		zp.setRemoveCif(false);
		ChemCompGroupFactory.setChemCompProvider(zp);

		long startTime = System.currentTimeMillis();
		StructureIO.getStructure("4HHM");
		long finishTime = System.currentTimeMillis();
		s_logger.info("ZipChemComp time: "+(finishTime-startTime)+ " ms");

		// Not wanted here for testing, but useful for cleaning up downloaded .cif.gz files.
		// ZipChemCompProvider.purgeTempFiles(pdbdir.toString());
	}
	
	@Test
	public void testGetOneLetterCode() throws Exception {
		
		String oneLetter;

		oneLetter = ChemCompGroupFactory.getOneLetterCode(ChemCompGroupFactory.getChemComp("ALA"));
		assertEquals("A", oneLetter);

		// single parent, we should return parent
		oneLetter = ChemCompGroupFactory.getOneLetterCode(ChemCompGroupFactory.getChemComp("MSE"));
		assertEquals("M", oneLetter);

		// multiparent case, we should return ?
		oneLetter = ChemCompGroupFactory.getOneLetterCode(ChemCompGroupFactory.getChemComp("OIM"));
		assertEquals("?", oneLetter);

	}

}
