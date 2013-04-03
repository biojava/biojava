/*
 *                    PDB web development code
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
 *
 * Created on Aug 5, 2009
 * Created by ap3
 *
 */

package org.biojava.bio.structure.secstruc;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.HetatomImpl;

public class SecStrucGroup
extends HetatomImpl
implements Group {

    /**
    *
    */
   private static final long serialVersionUID = 313490286720467714L;
    Atom N;
    Atom CA;
    Atom C;
    Atom O;
    Atom H;
    Group original;

    public SecStrucGroup(){
        super();
    }

    public String toString(){

        StringBuffer str = new StringBuffer("SecStrucGroup ");
        str.append(residueNumber);
        str.append(" ");
        str.append(pdb_name);
        str.append(" ");
        str.append(pdb_flag);
        if (pdb_flag) {
            str.append(" atoms: ");
            str.append(atoms.size());
        }


        return str.toString() ;

    }


    public Group getOriginal()
   {
      return original;
   }

   public void setOriginal(Group original)
   {
      this.original = original;
   }

   public Atom getC() {
        return C;
    }
    public void setC(Atom c) {
        addAtom(c);
        C = c;
    }
    public Atom getCA() {
        return CA;
    }
    public void setCA(Atom ca) {
        addAtom(ca);
        CA = ca;
    }
    public Atom getH() {
        return H;
    }
    public void setH(Atom h) {
        addAtom(h);
        H = h;
    }
    public Atom getN() {
        return N;
    }
    public void setN(Atom n) {
        addAtom(n);
        N = n;
    }
    public Atom getO() {
        return O;
    }
    public void setO(Atom o) {
        addAtom(o);
        O = o;
    }



}