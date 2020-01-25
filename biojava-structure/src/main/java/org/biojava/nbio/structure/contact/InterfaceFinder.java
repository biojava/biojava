package org.biojava.nbio.structure.contact;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.xtal.CrystalTransform;
import org.biojava.nbio.structure.xtal.SpaceGroup;

import javax.vecmath.Point3d;
import java.util.List;

/**
 * A class containing methods to find interfaces in a given structure.
 * @author Jose Duarte
 * @since 5.4.0
 */
public class InterfaceFinder {

    public static final double DEFAULT_CONTACT_CUTOFF = 6;

    private static final CrystalTransform IDENTITY_TRANSFORM = new CrystalTransform((SpaceGroup) null);
    private static final boolean INCLUDE_HETATOMS = true;

    private Structure structure;
    private double cutoff;

    private BoundingBox[] boundingBoxes;

    public InterfaceFinder(Structure structure) {
        this.structure = structure;
        this.cutoff = DEFAULT_CONTACT_CUTOFF;
    }

    /**
     * Set the contact distance cutoff.
     * @param cutoff the distance value in Angstroms
     */
    public void setCutoff(double cutoff) {
        this.cutoff = cutoff;
    }

    /**
     * Find all inter polymer-chain interfaces in the structure.
     * Two chains will be considered in contact if at least a pair of atoms (one from each chain) is within the
     * contact cutoff.
     * @return the list of all interfaces
     */
    public StructureInterfaceList getAllInterfaces() {
        initBoundingBoxes();

        StructureInterfaceList list = new StructureInterfaceList();

        List<Chain> polyChains = structure.getPolyChains();
        for (int i = 0; i<polyChains.size(); i++) {
            for (int j = i + 1; j<polyChains.size(); j++) {
                if (! boundingBoxes[i].overlaps(boundingBoxes[j], cutoff)) {
                    continue;
                }
                StructureInterface interf = calcInterface(polyChains.get(i), polyChains.get(j));
                if (interf!=null) {
                    list.add(interf);
                }
            }
        }
        return list;
    }

    private void initBoundingBoxes() {
        List<Chain> polyChains = structure.getPolyChains();
        boundingBoxes = new BoundingBox[polyChains.size()];
        for (int i = 0; i<polyChains.size(); i++) {
            Atom[] atoms = StructureTools.getAllNonHAtomArray(polyChains.get(i), INCLUDE_HETATOMS);
            Point3d[] points = Calc.atomsToPoints(atoms);
            BoundingBox bb = new BoundingBox(points);
            boundingBoxes[i] = bb;
        }
    }

    private StructureInterface calcInterface(Chain chain1, Chain chain2) {
        AtomContactSet graph = StructureTools.getAtomsInContact(chain1, chain2, cutoff, INCLUDE_HETATOMS);

        StructureInterface interf = null;
        if (graph.size()>0) {
            interf = new StructureInterface(
                    StructureTools.getAllNonHAtomArray(chain1, INCLUDE_HETATOMS), StructureTools.getAllNonHAtomArray(chain2, INCLUDE_HETATOMS),
                    chain1.getName(), chain2.getName(),
                    graph,
                    IDENTITY_TRANSFORM, IDENTITY_TRANSFORM);
        }

        return interf;
    }
}
