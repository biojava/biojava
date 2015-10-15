package org.biojava.nbio.structure.secstruc;

/**
 * Container that represents a beta Bridge between two residues.
 * It contains the two partner indices and the type of the bridge.
 * For consistency, partner1 is always the small index.
 * 
 * @author Aleix Lafita
 *
 */
public class BetaBridge {
	
    BridgeType type;
    int partner1;
    int partner2;
    
    public BetaBridge(int i, int j, BridgeType t) {
		partner1 = Math.min(i, j);
		partner2 = Math.max(i, j);
		type = t;
	}
    
    @Override
    public boolean equals(Object o){
    	
    	if (!(o instanceof BetaBridge)) return false;
    	
    	BetaBridge b = (BetaBridge) o;
    	if (type != b.type) return false;
    	if (partner1 != b.partner1) return false;
    	if (partner2 != b.partner2) return false;
    	return true;
    }

	public BridgeType getType() {
		return type;
	}

	public int getPartner1() {
		return partner1;
	}

	public int getPartner2() {
		return partner2;
	}
    
}
