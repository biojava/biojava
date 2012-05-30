package demo;


import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;

import org.biojava3.structure.StructureIO;

public class DemoShowBiolAssembly {

	public static void main(String[] args){

		

		try{
			// see also: http://www.pdb.org/pdb/101/static101.do?p=education_discussion/Looking-at-Structures/bioassembly_tutorial.html
			// good examples: 1stp 1gav 1hv4 1hho 7dfr 3fad  1qqp
			
			Structure bioAssembly = StructureIO.getBiologicalAssembly("1stp",0);
			System.out.println(bioAssembly);
			
			StructureAlignmentJmol jmolPanel = new StructureAlignmentJmol();
			//jmolPanel.evalString("set autobond=false");
			jmolPanel.setStructure(bioAssembly);

			// send some commands to Jmol
			jmolPanel.evalString("select * ; color structure ; spacefill off; wireframe off; backbone off; cartoon on; select ligands ; spacefill 0.4; color cpk;");
			
			System.out.println("done!");
			//quaternaryBuilder.getBioUnitTransformationList()


		} catch (Exception e) {
			e.printStackTrace();
		}	
	}
}
