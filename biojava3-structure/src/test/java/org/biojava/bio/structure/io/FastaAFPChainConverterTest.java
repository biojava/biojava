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
	public void testCp() throws IOException, StructureException {
		
// initial unaligned:
//alfdynatgdtefdspakqgwmqdntnngsgvltnadgmpawlvqgiggraqwtyslstnqhaqassfgwrmttemkvlsggmitnyyangtqrvlpiisldssgnlvvefegqtgrtvlatgtaateyhkfelvflpgsnpsasfyfdgklirdniqptaskQNMIVWGNGSS
//--------------------------------------------------------------------------------------------kirivdgaanqiqvadgsrkyvvtlsidesgglvanlngvsapiilqsehakvhsfhdyelqysalnhtttLFVDGQQITTW
		// 753 - 393 = 360
// total length: 1094
// 1094 - 393 = 701 = actual cp site
// have: 892
//alfdynatgdtefdspakqgwmqdntnngsgvltnadgmpawlvqgiggraqwtyslstnqhaqassfgwrmttemkvlsggmitnyyangtqrvlpiisldssgnlvvefegqtgrtvlatgtaateyhkfelvflpgsnpsasfyfdgklirdniqptaskQNMIVWGNGSSntdgvaayrdikfei------------------------------------------------------------------------------------------------------------------QGDVIf------------RGPDRIPSIVASsvTPGVVTAFAEKRVGGgdpgalsntNDIITRTSRDGGITWDTELNLTEQinvsdeFDFSDPRPIYDPs---SNTVLVSYARWPtdaaqngdrikpwmpNGIFYSVYDVASgnWQAPIDVTdqvkersfqiagwggselyrrntslnsqqdwqsnakirivdgaanqiqvadgsrkyvvtlsidesgglvanlngvsapiilqsehakvhsfhdyelqysalnhtttlfvdgqqittwagevsqenniqfgnadaqidgrlhvqkivltqqghnlvefdafylaqqtpevekdleklgwtkiktgntmslygNASVNPGpgHGITLtrqqnisgsqNGRLIYPAIVLdrfFLNVMSIYSDDGgsnwq-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TGSTLpipfrwksssileTLEPSEADMVELQN--GDLLLTARLDFNQivngvny--SPRQQFLSKDGGITWSLLEANNANvfsnistgTVDASITRFEqsdgSHFLLFTNPQGnpagTNgr------------QNLGLWFSFDEG--VTWKGPIQ--LVNGasaysdiyqldsenaivivetdnsnmrilrmpitllkqklt
// so permuted by +360 should be:
//NDIITRTSRDGGITWDTELNLTEQinvsdeFDFSDPRPIYDPs---SNTVLVSYARWPtdaaqngdrikpwmpNGIFYSVYDVASgnWQAPIDVTdqvkersfqiagwggselyrrntslnsqqdwqsnakirivdgaanqiqvadgsrkyvvtlsidesgglvanlngvsapiilqsehakvhsfhdyelqysalnhtttlfvdgqqittwagevsqenniqfgnadaqidgrlhvqkivltqqghnlvefdafylaqqtpevekdleklgwtkiktgntmslygNASVNPGpgHGITLtrqqnisgsqNGRLIYPAIVLdrfFLNVMSIYSDDGgsnwq-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TGSTLpipfrwksssileTLEPSEADMVELQN--GDLLLTARLDFNQivngvny--SPRQQFLSKDGGITWSLLEANNANvfsnistgTVDASITRFEqsdgSHFLLFTNPQGnpagTNgr------------QNLGLWFSFDEG--VTWKGPIQ--LVNGasaysdiyqldsenaivivetdnsnmrilrmpitllkqklt
// but rotated back by -360 should be:
//DPRPIYDPs---SNTVLVSYARWPtdaaqngdrikpwmpNGIFYSVYDVASgnWQAPIDVTdqvkersfqiagwggselyrrntslnsqqdwqsnakirivdgaanqiqvadgsrkyvvtlsidesgglvanlngvsapiilqsehakvhsfhdyelqysalnhtttlfvdgqqittwagevsqenniqfgnadaqidgrlhvqkivltqqghnlvefdafylaqqtpevekdleklgwtkiktgntmslygNASVNPGpgHGITLtrqqnisgsqNGRLIYPAIVLdrfFLNVMSIYSDDGgsnwq-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TGSTLpipfrwksssileTLEPSEADMVELQN--GDLLLTARLDFNQivngvny--SPRQQFLSKDGGITWSLLEANNANvfsnistgTVDASITRFEqsdgSHFLLFTNPQGnpagTNgr------------QNLGLWFSFDEG--VTWKGPIQ--LVNGasaysdiyqldsenaivivetdnsnmrilrmpitllkqklt
		Structure structure = cache.getStructure("1w0p");
		//!  IS     AS    aas   nm  nma  angle  transl. idm nma/nr   Ts/nr   TM/nr Zc(Ts) Z1(Ts) Z1(Tm) S?   Ux       Uy       Uz
		// -393   42.4  334.3  249  249  180.6    0.21   0  0.3307  0.1357  0.2949  18.16  19.13  13.27  1  0.83108 -0.54728  0.09893 
		
		
		
/*
INPUT:
alfdynatgdtefdspakqgwmqdntnngsgvltnadgmpawlvqgiggraqwtyslstnqhaqassfgwrmttemkvlsggmitnyyangtqrvlpiisldssgnlvvefegqtgrtvlatgtaateyhkfelvflpgsnpsasfyfdgklirdniqptaskQNMIVWGNGSSntdgvaayrdikfei------------------------------------------------------------------------------------------------------------------QGDVIf------------RGPDRIPSIVASsvTPGVVTAFAEKRVGGgdpgalsntNDIITRTSRDGGITWDTELNLTEQinvsdeFDFSDPRPIYDPs---SNTVLVSYARWPtdaaqngdrikpwmpNGIFYSVYDVASgnWQAPIDVTdqvkersfqiagwggselyrrntslnsqqdwqsnakirivdgaanqiqvadgsrkyvvtlsidesgglvanlngvsapiilqsehakvhsfhdyelqysalnhtttlfvdgqqittwagevsqenniqfgnadaqidgrlhvqkivltqqghnlvefdafylaqqtpevekdleklgwtkiktgntmslygNASVNPGpgHGITLtrqqnisgsqNGRLIYPAIVLdrfFLNVMSIYSDDGgsnwq-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TGSTLpipfrwksssileTLEPSEADMVELQN--GDLLLTARLDFNQivngvny--SPRQQFLSKDGGITWSLLEANNANvfsnistgTVDASITRFEqsdgSHFLLFTNPQGnpagTNgr------------QNLGLWFSFDEG--VTWKGPIQ--LVNGasaysdiyqldsenaivivetdnsnmrilrmpitllkqklt
--------------------------------------------------------------------------------------------kirivdgaanqiqvadgsrkyvvtlsidesgglvanlngvsapiilqsehakvhsfhdyelqysalnhtttLFVDGQQITTWagevsqenniqfgnadaqidgrlhvqkivltqqghnlvefdafylaqqtpevekdleklgwtkiktgntmslygnasvnpgpghgitltrqqnisgsqngrliypaivldrfflnvmsiysddggsnwqTGSTLpipfrwksssileTLEPSEADMVEL--QNGDLLLTARLDFNQivngvny--SPRQQFLSKDGGITWSLLEANNANvfsnisTGTVDASITRFEqsdgSHFLLFTNPQGNpagtngr--------QNLGLWFSFDEG--VTWKGPIQlv---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------NGASAYS--DIYQLd---------SENAIVIVETD---NSNMRILRMPITllkqkltalfdynatgdtefdspakqgwmqdntnngsgvltnadgmpawlvqgiggraqwtyslstnqhaqassfgwrmttemkvlsggmitnyyangtqrvlpiisldssgnlvvefegqtgrtvlatgtaateyhkfelvflpgsnpsasfyfdgklirdniqptaskqnmivwgngssntdgvaayrdikfeiQGDVIf------------RGPDRIPSIVASSVtpGVVTAFAEKRVGGgdpgalsntNDIITRTSRDGGITWDTELNLTEQinvsdefdFSDPRPIYDPs---SNTVLVSYARW----PTdaaqngdrikpwmpNGIFYSVYDVASgnWQAPIDVTdqVKERsfqiagwggselyrrntslnsqqdwqsna------------
OUTPUT:
===================================================================================================================================================================QNMIVWGNGSS=================================================================================================================================QGDVI=============RGPDRIPSIVAS==TPGVVTAFAEKRVGG=========NDIITRTSRDGGITWDTELNLTEQ======FDFSDPRPIYDP====SNTVLVSYARWP===============NGIFYSVYDVAS==WQAPIDVT===============================================================================================================================================================================================NASVNPG==HGITL==========NGRLIYPAIVL===FLNVMSIYSDDG====================================================================================================================================================================================================TGSTL=============TLEPSEADMVELQN==GDLLLTARLDFNQ=========SPRQQFLSKDGGITWSLLEANNAN========TVDASITRFE====SHFLLFTNPQG====TN==============QNLGLWFSFDEG==VTWKGPIQ==LVNG=========================================
===================================================================================================================================================================LFVDGQQITTW=================================================================================================================================TGSTL=============TLEPSEADMVEL==QNGDLLLTARLDFNQ=========SPRQQFLSKDGGITWSLLEANNAN======TGTV==================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================QGDVI=============RGPDRIPSIVASSV==GVVTAFAEKRVGG=========NDIITRTSRDGGITWDTELNLTEQ========FSDPRPIYDP====SNTVLVSYARW====PT==============NGIFYSVYDVAS==WQAPIDVT==VKER=========================================
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CP---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
*/
		
		
		
		// first aligned in permuted occurs at 918
		String first = ("alfdynatgdtefdspakqgwmqdntnngsgvltnadgmpawlvqgiggraqwtyslstnqhaqassfgwrmttemkvlsggmitnyyangtqrvlpiisldssgnlvvefegqtgrtvlatgtaateyhkfelvflpgsnpsasfyfdgklirdniqptaskQNMIVWGNGSSntdgvaayrdikfei------------------------------------------------------------------------------------------------------------------QGDVIf------------RGPDRIPSIVASsvTPGVVTAFAEKRVGGgdpgalsntNDIITRTSRDGGITWDTELNLTEQinvsdeFDFSDPRPIYDPs---SNTVLVSYARWPtdaaqngdrikpwmpNGIFYSVYDVASgnWQAPIDVTdqvkersfqiagwggselyrrntslnsqqdwqsnakirivdgaanqiqvadgsrkyvvtlsidesgglvanlngvsapiilqsehakvhsfhdyelqysalnhtttlfvdgqqittwagevsqenniqfgnadaqidgrlhvqkivltqqghnlvefdafylaqqtpevekdleklgwtkiktgntmslygNASVNPGpgHGITLtrqqnisgsqNGRLIYPAIVLdrfFLNVMSIYSDDGgsnwq-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TGSTLpipfrwksssileTLEPSEADMVELQN--GDLLLTARLDFNQivngvny--SPRQQFLSKDGGITWSLLEANNANvfsnistgTVDASITRFEqsdgSHFLLFTNPQGnpagTNgr------------QNLGLWFSFDEG--VTWKGPIQ--LVNGasaysdiyqldsenaivivetdnsnmrilrmpitllkqklt");
		String second =   ("--------------------------------------------------------------------------------------------kirivdgaanqiqvadgsrkyvvtlsidesgglvanlngvsapiilqsehakvhsfhdyelqysalnhtttLFVDGQQITTWagevsqenniqfgnadaqidgrlhvqkivltqqghnlvefdafylaqqtpevekdleklgwtkiktgntmslygnasvnpgpghgitltrqqnisgsqngrliypaivldrfflnvmsiysddggsnwqTGSTLpipfrwksssileTLEPSEADMVEL--QNGDLLLTARLDFNQivngvny--SPRQQFLSKDGGITWSLLEANNANvfsnisTGTVDASITRFEqsdgSHFLLFTNPQGNpagtngr--------QNLGLWFSFDEG--VTWKGPIQlv---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------NGASAYS--DIYQLd---------SENAIVIVETD---NSNMRILRMPITllkqkltalfdynatgdtefdspakqgwmqdntnngsgvltnadgmpawlvqgiggraqwtyslstnqhaqassfgwrmttemkvlsggmitnyyangtqrvlpiisldssgnlvvefegqtgrtvlatgtaateyhkfelvflpgsnpsasfyfdgklirdniqptaskqnmivwgngssntdgvaayrdikfeiQGDVIf------------RGPDRIPSIVASSVtpGVVTAFAEKRVGGgdpgalsntNDIITRTSRDGGITWDTELNLTEQinvsdefdFSDPRPIYDPs---SNTVLVSYARW----PTdaaqngdrikpwmpNGIFYSVYDVASgnWQAPIDVTdqVKERsfqiagwggselyrrntslnsqqdwqsna------------");

		int y = 0;
		for (int i = 0; i < second.length(); i++) if (second.charAt(i) == '-' || Character.isLowerCase(second.charAt(i))) y++;
		System.out.println("y=" + y);
		
		AFPChain afpChain = FastaAFPChainConverter.cpFastaToAfpChain(first, second, structure, -393);
		String xml = AFPChainXMLConverter.toXML(afpChain);
		System.out.println(xml);
	}
	
//	@Test
	public void testIncomplete() throws IOException, StructureException {
		Structure s1 = cache.getStructure("1w0p");
		Structure s2 = cache.getStructure("1qdm");
		ProteinSequence seq1 = new ProteinSequence("GWGG----SEL--YRRNTSLNS--QQDW-------QSNAKIRIVDGAA-----NQIQ");
		ProteinSequence seq2 = new ProteinSequence("WMQNQLAQNKT--QDLILDYVNQLCNRL---PSPMESAV----DCGSLGSMPDIEFT");
		AFPChain afpChain = FastaAFPChainConverter.fastaToAfpChain(seq1, seq2, s1, s2);
		assertEquals("Wrong number of EQRs", 33, afpChain.getNrEQR());
		String xml = AFPChainXMLConverter.toXML(afpChain);
		System.out.println(xml);
		File expected = new File("src/test/resources/1w0p_1qdm.xml");
		File x = File.createTempFile("1w0p_1qdm_output", "xml.tmp");
		x.deleteOnExit();
		BufferedWriter bw = new BufferedWriter(new FileWriter(x));
		bw.write(xml);
		bw.close();
		assertTrue("AFPChain is wrong", compareXml(expected, x));
	}

}
