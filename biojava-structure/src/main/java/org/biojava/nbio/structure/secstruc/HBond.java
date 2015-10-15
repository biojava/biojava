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
package org.biojava.nbio.structure.secstruc;

/**
 * Container that represents a hidrogen bond. 
 * It contains the energy of the bond in cal/mol and the partner index.
 * 
 * @author Andreas Prlic
 * @author Aleix Lafita
 *
 */
public class HBond{
	
    private double energy;
    private int partner;
    
    public HBond() {
        energy = 0;
        partner = 0;
    }
    
    public HBond(HBond o){
    	this.energy = o.energy;
    	this.partner = o.partner;
    }
    
    @Override
	public HBond clone(){
        return new HBond(this);
    }
    
    @Override
	public String toString(){
        return partner+" | "+(energy/1000.0);
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
