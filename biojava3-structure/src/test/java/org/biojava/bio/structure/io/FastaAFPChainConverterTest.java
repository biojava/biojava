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
 * Created on 2013-05-28
 * Created by Douglas Myers-Turnbull
 *
 * @since 3.0.6
 */
package org.biojava.bio.structure.io;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFP;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava3.core.sequence.ProteinSequence;
import org.custommonkey.xmlunit.DetailedDiff;
import org.custommonkey.xmlunit.Diff;
import org.custommonkey.xmlunit.Difference;
import org.custommonkey.xmlunit.XMLUnit;
import org.custommonkey.xmlunit.examples.RecursiveElementNameAndTextQualifier;
import org.junit.Before;
import org.junit.Test;
import org.xml.sax.SAXException;


/**
 * A test for {@link FastaAFPChainConverter}.
 * @author dmyersturnbull
 *
 */
public class FastaAFPChainConverterTest {

	static {
		XMLUnit.setIgnoreWhitespace(true);
		XMLUnit.setIgnoreComments(true);
		XMLUnit.setIgnoreAttributeOrder(true);
	}

	public static void printDetailedDiff(Diff diff, PrintStream ps) {
		DetailedDiff detDiff = new DetailedDiff(diff);
		for (Object object : detDiff.getAllDifferences()) {
			Difference difference = (Difference) object;
			ps.println(difference);
		}
	}

	/**
	 * Compares two XML files without regard to the order of elements or attributes, and ignoring any element named \"releaseDate\".
	 * @return Whether the files are \"similar\"
	 */
	public static boolean compareXml(File expectedFile, File actualFile) {
		try {
			FileReader expectedFr = new FileReader(expectedFile);
			FileReader actualFr = new FileReader(actualFile);
			Diff diff = new Diff(expectedFr, actualFr);
			// ignore order
			// look at element, id, and weight (weight is a nested element)
			diff.overrideElementQualifier(new RecursiveElementNameAndTextQualifier());
			final boolean isSimilar = diff.similar();
			if (!isSimilar) printDetailedDiff(diff, System.err);
			expectedFr.close();
			actualFr.close();
			return isSimilar;
		} catch (IOException e) {
			throw new RuntimeException(e);
		} catch (SAXException e) {
			throw new RuntimeException(e);
		}
	}

	private AtomCache cache;
	
	@Before
	public void setUp() {
		cache = new AtomCache();
	}
	
	@Test
	public void testCpAsymmetric() throws IOException, StructureException {
		Structure structure = cache.getStructure("1w0p");
		String first = ("alfdynatgdtefdspakqgwmqdntnngsgvltnadgmpawlvqgiggraqwtyslstnqhaqassfgwrmttemkvlsggmitnyyangtqrvlpiisldssgnlvvefegqtgrtvlatgtaateyhkfelvflpgsnpsasfyfdgklirdniqptaskQNMIVWGNGSSntdgvaayrdikfei------------------------------------------------------------------------------------------------------------------QGDVIf------------RGPDRIPSIVASsvTPGVVTAFAEKRVGGgdpgalsntNDIITRTSRDGGITWDTELNLTEQinvsdeFDFSDPRPIYDPs---SNTVLVSYARWPtdaaqngdrikpwmpNGIFYSVYDVASgnWQAPIDVTdqvkersfqiagwggselyrrntslnsqqdwqsnakirivdgaanqiqvadgsrkyvvtlsidesgglvanlngvsapiilqsehakvhsfhdyelqysalnhtttlfvdgqqittwagevsqenniqfgnadaqidgrlhvqkivltqqghnlvefdafylaqqtpevekdleklgwtkiktgntmslygNASVNPGpgHGITLtrqqnisgsqNGRLIYPAIVLdrfFLNVMSIYSDDGgsnwq-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TGSTLpipfrwksssileTLEPSEADMVELQN--GDLLLTARLDFNQivngvny--SPRQQFLSKDGGITWSLLEANNANvfsnistgTVDASITRFEqsdgSHFLLFTNPQGnpagTNgr------------QNLGLWFSFDEG--VTWKGPIQ--LVNGasaysdiyqldsenaivivetdnsnmrilrmpitllkqklt");
		String second =   ("--------------------------------------------------------------------------------------------kirivdgaanqiqvadgsrkyvvtlsidesgglvanlngvsapiilqsehakvhsfhdyelqysalnhtttLFVDGQQITTWagevsqenniqfgnadaqidgrlhvqkivltqqghnlvefdafylaqqtpevekdleklgwtkiktgntmslygnasvnpgpghgitltrqqnisgsqngrliypaivldrfflnvmsiysddggsnwqTGSTLpipfrwksssileTLEPSEADMVEL--QNGDLLLTARLDFNQivngvny--SPRQQFLSKDGGITWSLLEANNANvfsnisTGTVDASITRFEqsdgSHFLLFTNPQGNpagtngr--------QNLGLWFSFDEG--VTWKGPIQlv---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------NGASAYS--DIYQLd---------SENAIVIVETD---NSNMRILRMPITllkqkltalfdynatgdtefdspakqgwmqdntnngsgvltnadgmpawlvqgiggraqwtyslstnqhaqassfgwrmttemkvlsggmitnyyangtqrvlpiisldssgnlvvefegqtgrtvlatgtaateyhkfelvflpgsnpsasfyfdgklirdniqptaskqnmivwgngssntdgvaayrdikfeiQGDVIf------------RGPDRIPSIVASSVtpGVVTAFAEKRVGGgdpgalsntNDIITRTSRDGGITWDTELNLTEQinvsdefdFSDPRPIYDPs---SNTVLVSYARW----PTdaaqngdrikpwmpNGIFYSVYDVASgnWQAPIDVTdqVKERsfqiagwggselyrrntslnsqqdwqsna------------");
		AFPChain afpChain = FastaAFPChainConverter.cpFastaToAfpChain(first, second, structure, -393);
		assertEquals("Wrong TM-score", 0.2949, afpChain.getTMScore(), 0.001);
		assertEquals("Wrong RMSD", 3.605, afpChain.getTotalRmsdOpt(), 0.001);
	}

	@Test
	public void testCpSymmetric2() throws IOException,StructureException {
		String a = "--vRSLNCTLRDSQQ-KSLVMSG---PYELKALHLQgqdmeq-----QVVFSMSFVQGeesndkiPVALGLKEK-NLYLSSVLKdDKPTLQLESVdpknypkkkmekRFVFNKIEInn--KLEFESAQFpnWYISTSqAENmPVFLGGT----KGgqDITDFTMQFV---";
		String b = "esnDKIPVALGLKEKnLYLSSVLkddKPTLQLESVDpknypkkkmekRFVFNKIEINN-------KLEFESAQFpNWYISTSQA-ENMPVFLGGTkggqd-------ITDFTMQFVvrslNCTLRDSQQ--KSLVMS-GPY-ELKALHLqgqdME--QQVVFSMSFVqge";
		Structure structure = StructureTools.getStructure("31BI");
		AFPChain afpChain = FastaAFPChainConverter.cpFastaToAfpChain(a, b, structure, -101);
		assertEquals("Wrong TM-score", 0.6284, afpChain.getTMScore(), 0.001);
		assertEquals("Wrong RMSD", 2.50569, afpChain.getTotalRmsdOpt(), 0.001);
	}
	
	@Test
	public void testCpSymmetric1() throws IOException,StructureException {
		//cat 2GG6-best.fasta |tr -d \\n|pbcopy
		String a = "-SSRPATAR-KSSGLSGTVRIPGDKSISHRSFMFGGLA-SGETRITGLLEG-EDvINTGKAMQAMGARIRKEGd---------TWIIDGVgngglLAPEAPLD---FGNAATGCRLTMGLVGvydFDSTFIGDASLtkrp---MGRVLNPLREMGVQVKSEDgdrLPVTLRGPK---TPT---PITYRVpMASAQVKSAVLLAGLNTPGITTVIEpi---MTRDHTEKMLQGFGANLTVEtdadGVRTIRLEgRGKLTGQVIDVPGDPSSTAFPLVAALLVpGSDVTILNVLMNpTR-TGLILTLQEMGADIEVINprlaggedvaDLRVRSS-----TLKGVTVPedrAPSMIDEYPILAVAAAFAEGATVMNGLEELrvkesdrLSAVANGLKLNGVDCDEGE---TSLVVRGRPdgkGLGNasgAAVAT-HLDHRIAMSFLVMGLVSENPVTVDDatmIATSFPEFMDLMAGLGAKIELS---";
		String b = "dGVRTIRLEgRGKLTGQVIDVPGDPSSTAFPLVAALLVpGSDVTILNVLMNpTR-TGLILTLQEMGADIEVINprlaggedvaDLRVRSS-----TLKGVTVPedrAPSMIDEYPILAVAAAfaeGATVMNGLEELrvkesdrLSAVANGLKLNGVDCDEGE---TSLVVRGRPdgkGLGnasGAAVAT-HLDHRIAMSFLVMGLVSENPVTVDDatmiaTSFPEFMDLMAGLGAKIELS----SSRPATAR-KSSGLSGTVRIPGDKSISHRSFMFGGLA-SGETRITGLLEG-EDvINTGKAMQAMGARIRKEGd---------TWIIDGVgngglLAPEAPLD---FGNAATGCRLTMGLVGVYDFDSTFIGDASLtkrp---MGRVLNPLREMGVQVKSEDgdrLPVTLRGPK---TPTP---ITYRVpMASAQVKSAVLLAGLNTPGITTVIE---PIMTRDHTEKMLQGFGANLTVEtda";
		Structure structure = StructureTools.getStructure("2GG6");
		AFPChain afpChain = FastaAFPChainConverter.cpFastaToAfpChain(a, b, structure, -230); // 215
		assertEquals("Wrong TM-score", 0.7701, afpChain.getTMScore(), 0.001);
		assertEquals("Wrong RMSD", 3.035, afpChain.getTotalRmsdOpt(), 0.001);
	}
	
	@Test
	public void testFromFasta() throws IOException, StructureException {
		Structure s1 = cache.getStructure("1w0p");
		Structure s2 = cache.getStructure("1qdm");
		ProteinSequence seq1 = new ProteinSequence("GWGG----SEL--YRRNTSLNS--QQDW-------QSNAKIRIVDGAA-----NQIQ");
		ProteinSequence seq2 = new ProteinSequence("WMQNQLAQNKT--QDLILDYVNQLCNRL---PSPMESAV----DCGSLGSMPDIEFT");
		AFPChain afpChain = FastaAFPChainConverter.fastaToAfpChain(seq1, seq2, s1, s2);
		assertEquals("Wrong number of EQRs", 33, afpChain.getNrEQR());
		String xml = AFPChainXMLConverter.toXML(afpChain);
		File expected = new File("src/test/resources/1w0p_1qdm.xml");
		File x = File.createTempFile("1w0p_1qdm_output", "xml.tmp");
		x.deleteOnExit();
		BufferedWriter bw = new BufferedWriter(new FileWriter(x));
		bw.write(xml);
		bw.close();
		assertTrue("AFPChain is wrong", compareXml(expected, x));
	}

}
