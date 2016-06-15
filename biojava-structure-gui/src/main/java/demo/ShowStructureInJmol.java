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

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.StructureIO;


/** Demo how to load and display a structure in Jmol
 *
 * @author Andreas Prlic
 *
 */
public class ShowStructureInJmol {
	public static void main(String[] args){
		try {

			Structure struc = StructureIO.getStructure("1aoi");

			StructureAlignmentJmol jmolPanel = new StructureAlignmentJmol();

			jmolPanel.setStructure(struc);

			// send some RASMOL style commands to Jmol
			jmolPanel.evalString("select * ; color chain;");
			jmolPanel.evalString("select nucleic; cartoon on;");
			jmolPanel.evalString("select *; spacefill off; wireframe off; cartoon on;  ");
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
