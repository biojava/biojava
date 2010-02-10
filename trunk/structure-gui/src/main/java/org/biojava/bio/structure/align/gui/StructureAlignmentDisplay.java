package org.biojava.bio.structure.align.gui;

import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.jama.Matrix;

public class StructureAlignmentDisplay {

	public static StructureAlignmentJmol display(AFPChain afpChain, Atom[] ca1,
			Atom[] ca2, List<Group> hetatms, List<Group> nucs1,
			List<Group> hetatms2, List<Group> nucs2) throws StructureException {

		Group[] twistedGroups = new Group[ ca2.length];

		int i=-1;
		if ( afpChain.getBlockRotationMatrix().length > 0 ) {
			Matrix m =  afpChain.getBlockRotationMatrix()[0];
			Atom shift =  afpChain.getBlockShiftVector()[0];

			for (Atom a: ca2){
				i++;
				twistedGroups[i]=a.getParent();
				Calc.rotate(twistedGroups[i],m);
				Calc.shift(twistedGroups[i],shift);
			}
		}

		return DisplayAFP.display(afpChain, twistedGroups, ca1, ca2,hetatms, nucs1, hetatms2, nucs2);

	}
}
