package org.biojava.nbio.structure.basepairs;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;

import javax.vecmath.Matrix4d;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by luke on 7/22/17.
 */
public class TertiaryBasePairParameters extends BasePairParameters {

    public TertiaryBasePairParameters(Structure structure, boolean RNA, boolean removeDups) {
        super(structure, RNA, removeDups);
    }

    @Override
    public List<Group[]> findPairs(List<Chain> chains) {
        List<Group[]> result = new ArrayList<>();
        for (int i = 0; i < chains.size(); i++) {
            Chain c = chains.get(i);
            String sequence = c.getSeqResSequence();
            for (int j = 0; j < sequence.length(); j++) {
                for (int k = sequence.length()-1; k >= j + 4; k--) {
                    Group g1 = c.getSeqResGroup(j);
                    Group g2 = c.getSeqResGroup(k);
                    Integer type1 = map.get(g1.getPDBName());
                    Integer type2 = map.get(g2.getPDBName());
                    if (type1 == null || type2 == null) continue;
                    Atom a1 = g1.getAtom("C1'");
                    Atom a2 = g2.getAtom("C1'");
                    if (a1 == null || a2 == null) continue;
                    // C1'-C1' distance is one useful criteria
                    if (Math.abs(a1.getCoordsAsPoint3d().distance(a2.getCoordsAsPoint3d())-10.0) > 5.0) continue;
                    Group[] ga = new Group[] {g1, g2};
                    Matrix4d data = basePairReferenceFrame(ga);
                    // if the stagger is greater than 2 Å, it's not really paired.
                    if (Math.abs(pairParameters[5]) > 2.0) continue;
                    // if the propeller is ridiculous it's also not that good of a pair.
                    if (Math.abs(pairParameters[1]) > 60.0) {
                        continue;
                    }
                    result.add(ga);
                    pairingNames.add(useRNA ? baseListRNA[type1]+baseListRNA[type2]: baseListDNA[type1]+baseListDNA[type2]);
                }
            }
        }
        result.addAll(super.findPairs(chains));
        return result;
    }


}
