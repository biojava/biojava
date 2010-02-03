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
 * Created on 2010-01-21
 *
 */
package demo;

import org.biojava.bio.structure.align.gui.AlignmentGui;

/**
 * 
 * @author Andreas Prlic
 *
 */
public class AlignmentGuiDemo {

	//  /Users/ap3/WORK/PDB/hy/pdb2hyn.ent.gz
	//  /Users/ap3/WORK/PDB/zl/pdb1zll.ent.gz
	
	public static void main(String[] args) {
		System.setProperty("PDB_DIR","/tmp/");
	
		AlignmentGui.getInstance();
	}
}
