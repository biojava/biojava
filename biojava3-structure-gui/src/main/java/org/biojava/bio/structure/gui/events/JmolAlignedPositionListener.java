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
 * created at May 28, 2008
 */
package org.biojava.bio.structure.gui.events;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.StructurePairAligner;
import org.biojava.bio.structure.gui.BiojavaJmol;
import org.biojava.bio.structure.gui.util.AlignedPosition;



public class JmolAlignedPositionListener implements AlignmentPositionListener{

	BiojavaJmol parent;
	Atom[] ca1;
	Atom[] ca2;
	StructurePairAligner structurePairAligner;
	
	public JmolAlignedPositionListener(BiojavaJmol parent, StructurePairAligner alig){
		this.parent = parent;
		structurePairAligner = alig;
	}
	


	public void setStructure1(Structure structure1) {
		
		ca1 = structurePairAligner.getAlignmentAtoms(structure1);
	}

	

	public void setStructure2(Structure structure2) {
		
		ca2 = structurePairAligner.getAlignmentAtoms(structure2);
	}

	public void mouseOverPosition(AlignedPosition p) {
		
		//System.out.println("mouseoverposition " + p);
		
		int p1 = p.getPos1();
		int p2 = p.getPos2();
		String s = "select ";
		
		if ((p1 > ca1.length) || (p2 > ca2.length)){
			System.err.println("requsting atom out of bounds! " );
			return;
		}
		
		String pdbpos1 ="";
		String pdbpos2 = "";
		
		if ( p1 >-1) {
			Atom a = ca1[p1];
			Group parent = a.getGroup();
			Chain c = parent.getChain();
			pdbpos1 = parent.getResidueNumber().toString();
			//System.out.printlng"chainid 1 is >"+c.getName()+"<");
			if (! c.getChainID().equals( " ")) {
				pdbpos1 += ":" + c.getChainID();
			} 
				
			
			//System.out.println("1:" + parent);
			s += pdbpos1 +"/1";
		}
		
		if ( p2 >-1) {
			Atom a = ca2[p2];
			Group parent = a.getGroup();
			Chain c = parent.getChain();
			pdbpos2 = parent.getResidueNumber().toString();
			//System.out.println("2:" + parent);
			//System.out.println("chainid 2 is >"+c.getName()+"<");
			if (! c.getChainID().equals( " ")) {
				pdbpos2 += ":" + c.getChainID();
			}
			if ( p1 > -1)
				s +=",";
			s += pdbpos2 +"/2";
		}
		s+="; set display selected;";
		//System.out.println(s);
		parent.evalString(s);
		
		
	}

	public void positionSelected(AlignedPosition p) {
		mouseOverPosition(p);
		
	}

	public void rangeSelected(AlignedPosition start, AlignedPosition end) {
		
		
	}

	public void selectionLocked() {
		// TODO Auto-generated method stub
		
	}

	public void selectionUnlocked() {
		// TODO Auto-generated method stub
		
	}



	public void toggleSelection(AlignedPosition p) {
		// TODO Auto-generated method stub
		
	}

	
	
}
