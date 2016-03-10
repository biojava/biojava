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
package demo;
import org.biojava.nbio.aaproperties.IPeptideProperties;
import org.biojava.nbio.aaproperties.PeptidePropertiesImpl;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;

/**
 * Created by andreas on 8/9/14.
 */


public class BioJavaAADemo {


	public static void main(String[] args) throws CompoundNotFoundException {
		ProteinSequence pSequence = new ProteinSequence("VLSPADKTNVKAAWGKVGAHAG");

		IPeptideProperties pp = new PeptidePropertiesImpl();

		System.out.println("Peptide Properties: " + pp.getIsoelectricPoint(pSequence));
	}
}
