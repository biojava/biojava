package demo;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava3.structure.StructureIO;


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
