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


import org.biojava.nbio.structure.align.gui.AlignmentGui;
import org.biojava.nbio.structure.align.gui.MultipleAlignmentGUI;

/**
 * Get an instance of the two alignment GUIs.
 * Try to align 2hyn vs. 1zll, for example.
 *
 * @author Andreas Prlic
 *
 */
public class AlignmentGuiDemo {

	public static void main(String[] args) {

		AlignmentGui.getInstance();
		MultipleAlignmentGUI.getInstance();
	}
}
