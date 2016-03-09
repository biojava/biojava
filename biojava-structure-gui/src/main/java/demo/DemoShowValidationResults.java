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
 * created at Sep 18, 2013
 * Author: ap3
 */

package demo;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.validation.*;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import java.io.InputStream;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.List;
import java.util.zip.GZIPInputStream;

public class DemoShowValidationResults {

	public static void main(String[] args){
		//String pdbId ="3zjr";
		String pdbId ="3vtq";
		showPdbValidation(pdbId);
	}

	private static void showPdbValidation(String pdbId) {
		try {
			JAXBContext ctx = JAXBContext.newInstance(new Class[] {WwPDBValidationInformation.class});

			Unmarshaller um = ctx.createUnmarshaller();

			InputStream inStream = new GZIPInputStream(DemoShowValidationResults.class.getResourceAsStream("/"+pdbId+"-valdata.xml.gz"));

			WwPDBValidationInformation validationReport = (WwPDBValidationInformation) um.unmarshal(inStream);

			Entry entry = validationReport.getEntry();

			System.out.println(pdbId + " " + entry.getPDBRevisionNumber() +
					"\t Rfree: " + entry.getDCCRfree() +
					"\t Clashscore " + entry.getClashscore() +
					"\t % Ramachandran outliers: "  + entry.getPercentRamaOutliers() +
					"\t % RSRC outliers: " + entry.getPercentRSRZOutliers() );


			StructureAlignmentJmol jmolPanel = new StructureAlignmentJmol();

			Structure s = StructureIO.getStructure(pdbId);

			jmolPanel.setStructure(s);

			jmolPanel.evalString("select *; color grey ; cartoon off ; ");

			for (ModelledSubgroup subgroup:  validationReport.getModelledSubgroup()) {

				List<Clash> clashes = subgroup.getClash();

				String chainId = subgroup.getChain();
				//String resname = subgroup.getResname();
				String iCode = subgroup.getIcode();
				BigInteger resnum = subgroup.getResnum();
				//String altcode = subgroup.getAltcode();


				String pos = resnum.toString() ;
				if ( iCode !=null && iCode.length()>0 && (! iCode.equals(" ")))
					pos +="^" + iCode;
				pos +=":" + chainId;

				BigDecimal base = new BigDecimal(0.5);

				for (Clash clash : clashes){
					String clashatom = clash.getAtom();
					BigDecimal clashmag = clash.getClashmag();
					// pos1 icode A chain X should become:
					// 1^A:X
					// [MET]508:A.CA/1 #3918
					// insertion code: [ASP]1^A:A.CA/1 #2

					String clashj = pos + "." + clashatom;
					String jmols = " select " + clashj + "; color red; spacefill " + (base.add(clashmag)) + ";" ;
					System.out.println(jmols + " " + clashmag);
					jmolPanel.evalString(jmols);
				}


				for (AngleOutlier angleout : subgroup.getAngleOutlier()) {
					String atom0 = angleout.getAtom0();
					String atom1 = angleout.getAtom1();
					String atom2 = angleout.getAtom2();

					String anglej = "select " + pos + "." + atom0+"," +pos+"." + atom1 +"," + pos +"." + atom2+"; color wireframe blue; wireframe 0.5;";
					//System.out.println(anglej);
					jmolPanel.evalString(anglej);
				}

				for (BondOutlier bondout : subgroup.getBondOutlier()){
					String atom0 = bondout.getAtom0();
					String atom1 = bondout.getAtom1();
					String bondj = "select " + pos + "." + atom0+"," +pos+"." + atom1 +"; color wireframe green; wireframe 0.5;";
					jmolPanel.evalString(bondj);

				}
			}


		} catch (Exception e){
			e.printStackTrace();

		}

	}

}
