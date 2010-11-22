/*
 *                  BioJava development code
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
 * Created on Jan 7, 2006
 *
 */
package org.biojava.bio.structure.align.pairwise;


import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.jama.Matrix;


/** a pair of fragments of two protein structures
 * 
 * @author Andreas Prlic
 * @since 1.5
 * @version %I% %G%
 */
public class FragmentPair {

  
    int length;
    int pos1;
    int pos2;

    // parameter below may be used in different approaches
  
    int contacts;
    int cluster;
    double rms;
    int used;
    int covered;
    
    //filled if fragments are superimposed
    Matrix rot;
    Atom trans;

    //this unit vector indicates the rotation of j onto i
    Atom unitv;
    
    Atom center1;
    Atom center2;
    
    
    public FragmentPair(int length, int p1, int p2) {
        super();
        this.length = length ;
         pos1 = p1;
         pos2 = p2;
         
         contacts = 0;
         cluster = 0;
         rms = 0.0;
         used = 0;
         covered = 0;
         
         unitv = new AtomImpl();
         unitv.setX(0);
         unitv.setY(0);
         unitv.setZ(1);
         rot = null;
         trans = new AtomImpl();
         center1 = new AtomImpl();
         center2 = new AtomImpl();

    }
    public Object clone(){
        
        FragmentPair n = new FragmentPair(length,pos1,pos2);
        if ( center1 !=null)
            n.setCenter1((Atom)center1.clone());
        
        if ( center2 != null)
            n.setCenter2((Atom)center2.clone());
        
        n.setCluster(cluster);
        n.setContacts(contacts);
        n.setCovered(covered);
        n.setRms(rms);
        n.setLength(length);
        n.setRot((Matrix)rot.clone());
        n.setUnitv((Atom)unitv.clone());
        
        return n;
    }
    public int getCluster() {
        return cluster;
    }

    public void setCluster(int cluster) {
        this.cluster = cluster;
    }

    public int getContacts() {
        return contacts;
    }

    public void setContacts(int contacts) {
        this.contacts = contacts;
    }

    public int getCovered() {
        return covered;
    }

    public void setCovered(int covered) {
        this.covered = covered;
    }

    public int getLength() {
        return length;
    }

    public void setLength(int length) {
        this.length = length;
    }

    public int getPos1() {
        return pos1;
    }

    public void setPos1(int pos1) {
        this.pos1 = pos1;
    }

    public int getPos2() {
        return pos2;
    }

    public void setPos2(int pos2) {
        this.pos2 = pos2;
    }

    public double getRms() {
        return rms;
    }

    public void setRms(double rms) {
        this.rms = rms;
    }

    public Matrix getRot() {
        return rot;
    }

    public void setRot(Matrix rot) {
        this.rot = rot;
    }

    public Atom getTrans() {
        return trans;
    }

    public void setTrans(Atom trans) {
        this.trans = trans;
    }

    public Atom getUnitv() {
        return unitv;
    }

    public void setUnitv(Atom unitv) {
        this.unitv = unitv;
    }

    public int getUsed() {
        return used;
    }

    public void setUsed(int used) {
        this.used = used;
    }

    public Atom getCenter1() {
        return center1;
    }

    public void setCenter1(Atom center1) {
        this.center1 = center1;
    }

    public Atom getCenter2() {
        return center2;
    }

    public void setCenter2(Atom center2) {
        this.center2 = center2;
    }
    
    public String toString() {
    	return String.format("Fragment (%d,%d) len %d", pos1, pos2, length);
    }
}
