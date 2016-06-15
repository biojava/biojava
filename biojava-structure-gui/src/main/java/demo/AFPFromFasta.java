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
package demo;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.AfpChainWriter;
import org.biojava.nbio.structure.align.util.AFPChainScorer;
import org.biojava.nbio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.nbio.structure.io.FastaAFPChainConverter;
import org.biojava.nbio.structure.io.FastaStructureParser;

import java.io.IOException;

/**
 * Demo displaying a structural alignment from a FASTA file using {@link FastaAFPChainConverter}.
 *
 * @author dmyerstu
 * @see {@link DemoAlignmentFromFasta} Also demonstrates the display of {@link StructureAlignment StructureAlignments} from FASTA sequences, but does so using the more general
 *      {@link FastaStructureParser}
 */
public class AFPFromFasta {

	public static void main(String[] args) throws IOException, StructureException, Exception {
		Structure structure1 = StructureTools.getStructure("1w0p");
		Structure structure2 = StructureTools.getStructure("1w0p");
		String first = "alfdynatgdtefdspakqgwmqdntnngsgvltnadgmpawlvqgiggraqwtyslstnqhaqassfgwrmttemkvlsggmitnyyangtqrvlpiisldssgnlvvefegqtgrtvlatgtaateyhkfelvflpgsnpsasfyfdgklirdniqptaskQNMIVWGNGSSntdgvaayrdikfei------------------------------------------------------------------------------------------------------------------QGDVIf------------RGPDRIPSIVASsvTPGVVTAFAEKRVGGgdpgalsntNDIITRTSRDGGITWDTELNLTEQinvsdeFDFSDPRPIYDPs---SNTVLVSYARWPtdaaqngdrikpwmpNGIFYSVYDVASgnWQAPIDVTdqvkersfqiagwggselyrrntslnsqqdwqsnakirivdgaanqiqvadgsrkyvvtlsidesgglvanlngvsapiilqsehakvhsfhdyelqysalnhtttlfvdgqqittwagevsqenniqfgnadaqidgrlhvqkivltqqghnlvefdafylaqqtpevekdleklgwtkiktgntmslygNASVNPGpgHGITLtrqqnisgsqNGRLIYPAIVLdrfFLNVMSIYSDDGgsnwq-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TGSTLpipfrwksssileTLEPSEADMVELQN--GDLLLTARLDFNQivngvny--SPRQQFLSKDGGITWSLLEANNANvfsnistgTVDASITRFEqsdgSHFLLFTNPQGnpagTNgr------------QNLGLWFSFDEG--VTWKGPIQ--LVNGasaysdiyqldsenaivivetdnsnmrilrmpitllkqklt";
		String second = "--------------------------------------------------------------------------------------------kirivdgaanqiqvadgsrkyvvtlsidesgglvanlngvsapiilqsehakvhsfhdyelqysalnhtttLFVDGQQITTWagevsqenniqfgnadaqidgrlhvqkivltqqghnlvefdafylaqqtpevekdleklgwtkiktgntmslygnasvnpgpghgitltrqqnisgsqngrliypaivldrfflnvmsiysddggsnwqTGSTLpipfrwksssileTLEPSEADMVEL--QNGDLLLTARLDFNQivngvny--SPRQQFLSKDGGITWSLLEANNANvfsnisTGTVDASITRFEqsdgSHFLLFTNPQGNpagtngr--------QNLGLWFSFDEG--VTWKGPIQlv---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------NGASAYS--DIYQLd---------SENAIVIVETD---NSNMRILRMPITllkqkltalfdynatgdtefdspakqgwmqdntnngsgvltnadgmpawlvqgiggraqwtyslstnqhaqassfgwrmttemkvlsggmitnyyangtqrvlpiisldssgnlvvefegqtgrtvlatgtaateyhkfelvflpgsnpsasfyfdgklirdniqptaskqnmivwgngssntdgvaayrdikfeiQGDVIf------------RGPDRIPSIVASSVtpGVVTAFAEKRVGGgdpgalsntNDIITRTSRDGGITWDTELNLTEQinvsdefdFSDPRPIYDPs---SNTVLVSYARW----PTdaaqngdrikpwmpNGIFYSVYDVASgnWQAPIDVTdqVKERsfqiagwggselyrrntslnsqqdwqsna------------";
		AFPChain afpChain = FastaAFPChainConverter.cpFastaToAfpChain(first, second, structure1, -393);
		Atom[] ca1 = StructureTools.getAtomCAArray(structure1);
		Atom[] ca2 = StructureTools.getAtomCAArray(structure2);
		String xml = AFPChainXMLConverter.toXML(afpChain);
		System.out.println(xml);
		double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
		afpChain.setTMScore(tmScore);
		System.out.println(AfpChainWriter.toScoresList(afpChain));
		StructureAlignmentDisplay.display(afpChain, ca1, ca2);
	}

}
