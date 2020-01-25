package org.biojava.nbio.structure.contact;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.xtal.CrystalTransform;
import org.biojava.nbio.structure.xtal.SpaceGroup;

import java.util.List;

public class InterfaceFinder {

    public static final double DEFAULT_CONTACT_CUTOFF = 6;

    private static final CrystalTransform IDENTITY_TRANSFORM = new CrystalTransform(SpaceGroup.parseSpaceGroup("P1"));
    private static final boolean INCLUDE_HETATOMS = true;

    private Structure structure;
    private double cutoff;

    public InterfaceFinder(Structure structure) {
        this.structure = structure;
        this.cutoff = DEFAULT_CONTACT_CUTOFF;
    }

    public void setCutoff(double cutoff) {
        this.cutoff = cutoff;
    }

    public StructureInterfaceList getAllInterfaces() {
        StructureInterfaceList list = new StructureInterfaceList();

        List<Chain> polyChains = structure.getPolyChains();
        for (int i = 0; i<polyChains.size(); i++) {
            for (int j = i + 1; j<polyChains.size(); j++) {
                StructureInterface interf = calcInterface(polyChains.get(i), polyChains.get(j));
                if (interf!=null) {
                    list.add(interf);
                }
            }
        }
        return list;
    }

    private StructureInterface calcInterface(Chain chain1, Chain chain2) {
        AtomContactSet graph = StructureTools.getAtomsInContact(chain1, chain2, cutoff, INCLUDE_HETATOMS);

        StructureInterface interf = null;
        if (graph.size()>0) {
            interf = new StructureInterface(
                    StructureTools.getAllAtomArray(chain1), StructureTools.getAllAtomArray(chain2),
                    chain1.getName(), chain2.getName(),
                    graph,
                    IDENTITY_TRANSFORM, IDENTITY_TRANSFORM);
        }

        return interf;
    }
}
