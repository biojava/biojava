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


public class HBond{
    double energy;
    int partner;
    
    public HBond() {
        energy = 0;
        partner = 0;
    }
    
    public Object clone(){
        
        HBond n = new HBond();
        n.setEnergy(energy);
        n.setPartner(partner);
        
        return n;
        
    }
    
    public String toString(){
        StringBuffer buf = new StringBuffer();
        
        buf.append(partner+"|"+(energy/1000.0));
        return buf.toString();
            
        
    }
    public double getEnergy() {
        return energy;
    }
    public void setEnergy(double energy) {
        this.energy = energy;
    }

    public int getPartner() {
        return partner;
    }

    public void setPartner(int partner) {
        this.partner = partner;
    }
    
}
